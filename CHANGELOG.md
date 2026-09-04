# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
