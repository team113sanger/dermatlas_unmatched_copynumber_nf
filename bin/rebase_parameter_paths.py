#!/usr/bin/env python3
"""Rewrite every file path in a fur-cnvkit parameters.json to a bare basename.

`fur_cnvkit cnvkit_static_files` records the absolute path of every input it was
given - the BAMs, the metadata manifest, the reference FASTA, the refFlat file and
the four static BED/gene files. Under Nextflow those absolute paths point into the
task directory of the process that wrote the file, which no longer exists by the time
a downstream process reads it.

fur-cnvkit resolves the recorded paths with `Path(...)`, so a bare filename resolves
against the reading process's working directory. Rewriting each path to its basename
therefore makes the parameter file portable between tasks, provided every process that
consumes it also stages the files it needs (see the `path` inputs of CNVKIT_CALL_COPY_NUMBER
and CNVKIT_POOLED_REFERENCE).

Only the keys listed in PATH_KEYS are touched, so `unplaced_contig_prefixes` and
`metadata_columns` - which hold contig prefixes and column names, not paths - are left
exactly as fur-cnvkit wrote them.
"""

import argparse
import json
import os
import sys
import typing as t

# Scalar path values, and list-of-path values, as written by
# fur_cnvkit.generate_cnvkit_static_files.generate_parameter_file.
PATH_KEYS: t.Tuple[str, ...] = (
    "reference_fasta",
    "baitset_bed",
    "refflat_file",
    "sample_metadata_xlsx",
    "sample_metadata_file",
    "access_bed",
    "targets_bed",
    "antitargets_bed",
    "baitset_genes_file",
)
PATH_LIST_KEYS: t.Tuple[str, ...] = ("all_bams", "tumour_bams", "normal_bams")


def rebase(parameters: t.Dict[str, t.Any]) -> t.Dict[str, t.Any]:
    for key in PATH_KEYS:
        if key in parameters and parameters[key] is not None:
            parameters[key] = os.path.basename(str(parameters[key]))
    for key in PATH_LIST_KEYS:
        if key in parameters and parameters[key] is not None:
            parameters[key] = [os.path.basename(str(path)) for path in parameters[key]]
    return parameters


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "-i", required=True, help="parameters.json to read")
    parser.add_argument("--output", "-o", required=True, help="parameters.json to write")
    args = parser.parse_args()

    with open(args.input) as handle:
        parameters = json.load(handle)

    # A parameter file with no BAMs means the metadata manifest matched none of the
    # input BAMs - usually a sample-id or tumour/normal column mismatch. Every
    # downstream step would then run to completion and produce nothing, so fail here.
    if not parameters.get("all_bams"):
        print(
            f"ERROR: {args.input} lists no BAM files. Check that the sample-id column "
            "of the metadata manifest matches the BAM filenames.",
            file=sys.stderr,
        )
        return 1

    with open(args.output, "w") as handle:
        json.dump(rebase(parameters), handle, indent=4)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
