# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Keywords

As of version 0.2.0 the following *keywords* are used at the start of each
changelog entry to indicate the impact of the change:

- **REPRODUCIBILITY** - a change to the pipeline's scientific processing that
  may cause the same input data to produce different scientific outputs or
  results, including changes to algorithms, tolerances, randomisation,
  scientific functionality, or output formats.
- **ROBUSTNESS** - a fix or improvement to the pipeline's scientific
  functionality that improves correctness, reliability, or the range of inputs
  that can be processed, without intentionally changing the scientific results
  of an equivalent successful analysis.
- **INTEGRATION** - a change to how the pipeline integrates with other systems
  or infrastructure, without changing its scientific processing or results.

## [Unreleased]

### Added
- **INTEGRATION** - the launcher reports each successful run's work-dir usage to the
  Dermatlas website: `stats/resource-stats-<RUN_ID>.txt` now carries `pipeline_slug=`,
  and `on_pipeline_exit` feeds the file to `dermatlas-http cohort analysis-workdir-stats`
  (module-loaded via the new `DERMATLAS_HTTP_MODULE` setting, default `dermatlas-http`)
  after the stats are written and before the work dir is removed. Best effort: a failed
  report is a `NOTE:` with the command to replay, never a failed run. Requires
  dermatlas-web-client >= 0.6.2.

### Changed
- **INTEGRATION** - **Breaking:** `unmatched_copynumber.config` takes every location from
  the variable dermanager exports for it rather than rebuilding paths from a filename
  convention: `bam_path` from `BAMS_DIR`, `all_samples` and the `all_tumours` subcohort
  from `DNA_TUMOUR_LIST_ANALYSED_ALL`, the `one_tumour_per_patient` subcohort from
  `DNA_TUMOUR_LIST_ONE_TUMOUR_PER_PATIENT_ALL`, and `metadata_manifest` from
  `COHORT_METADATA_FILE` (a full path) instead of `${PROJECT_DIR}/metadata/${METADATA_FILE}`
  (a filename). `run_unmatched_copynumber.sh` checks all four at launch in place of
  `METADATA_FILE`, and the MANUAL ENVIRONMENT OVERRIDES block and the README's standalone
  contract list them: ten exports, not seven. A `source_me.sh` still exporting only
  `METADATA_FILE` now fails the launch with the missing variable named.
- **INTEGRATION** - `lib/Utils.groovy` and `.github/workflows/publish-assets.yml` are
  byte-identical to the `dermatlas_rnafusions_nf` 0.4.15 copies again (comment-only
  differences: pipeline-neutral docstrings, and a resolve-step comment that still said the
  rolling tag is force-moved).


## [0.1.1]

### Changed
- Updated pipeline slug for reporting via derm tracking


## [0.1.0]

Initial release. Ports the `unmatched_copynumber` wrapper scripts (0.4.0) to Nextflow.

### Added
- `main.nf` with two entry points: the default one calls a cohort against a pre-computed
  pooled normal reference; `-entry BUILD_REFERENCE` generates the static files and the
  pooled normal reference a baitset is called against.
- Per-subcohort analysis (`params.subcohorts`), replacing the single
  `GISTIC_INCLUDE_LIST` / `GISTIC_EXCLUDE_LIST` pair. Each subcohort gets its own
  filtered segments, GISTIC2 input `.seg`, penetrance plots and GISTIC2 run.
- Input validation before any task is scheduled: manifest columns, sample lists,
  subcohort definitions and the threshold string are all checked up front.
- Run reporting shared with the other Dermatlas pipelines - Slack notification and
  Dermatlas analysis-log record, both opt-in (`lib/Utils.groovy`).
- `assets/run_unmatched_copynumber.sh` launcher and `assets/unmatched_copynumber.config`,
  published as a `projectify` asset bundle on release.
- nf-test suite covering both entry points, the sample-exclusion path and the missing
  reference error.

### Changed
- The preliminary penetrance sweep runs one task per purity in parallel, rather than a
  serial loop over `purity_vs_log2ratio.tsv`.
- Re-calling `cn` runs one task per sample, and always runs `cnvkit call`: where the
  original script copied the file when the thresholds were unchanged, re-calling with the
  same thresholds reproduces it.
- `filter_cnvkit.py` and `visualise_penetrance.py` moved to `bin/`, so Nextflow puts them
  on `PATH` inside the fur-cnvkit container instead of them being invoked by absolute path
  through `singularity exec`.
- GISTIC2 input filenames are uniform across the two arms:
  `<prefix>_{filtered,unfiltered}_segments.seg` and
  `<prefix>_gistic_{filtered,unfiltered}_samples.txt`. The unfiltered arm previously used
  a different pattern from the filtered one.
- Samples excluded from a subcohort default to the hypersegmented samples this run
  identified, rather than to a path assembled from environment variables.
- The metadata manifest must be a TSV. The pipeline reads the sex column itself to
  re-call `cn` per sample, which it cannot do from an `.xlsx` sheet.
