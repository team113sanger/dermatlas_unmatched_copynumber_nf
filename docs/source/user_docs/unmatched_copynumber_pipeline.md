# Unmatched copy number calling

This page describes how unmatched (tumour-only) CNA calling is performed on WES cohorts in Dermatlas, and
what each step of `dermatlas_unmatched_copynumber_nf` does. It builds on
[fur-cnvkit](https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_cnvkit), a CNVkit extension developed in
the Adams lab at WTSI for the FUR project.

The pipeline replaces the numbered wrapper scripts of the `unmatched_copynumber` repository. If you have run
those, the mapping is:

| Was | Is now |
| --- | --- |
| `01_generate_static_files_wes_bed.sh` | `CNVKIT_STATIC_FILES` (only under `-entry BUILD_REFERENCE`) |
| `01_generate_parameters_files_wes_bed.sh` | `CNVKIT_PARAMETER_FILE` |
| `02_generate_pooled_ref_wes_bed.sh` | `CNVKIT_POOLED_REFERENCE` (only under `-entry BUILD_REFERENCE`) |
| `03_call_copy_numbers_wes_bed.sh` | `CNVKIT_CALL_COPY_NUMBER` |
| `04_generate_filtered_samples_list.sh` | `FIND_HYPERSEGMENTED_SAMPLES` |
| `05_generate_prelim_cn_plots.sh` | `PRELIMINARY_PENETRANCE_PLOT`, one task per purity |
| `06_re-call_cnvkit_cn.sh` | `CNVKIT_RECALL`, one task per sample |
| `07_filter_cns_segments.sh` | `FILTER_SEGMENTS` and `PENETRANCE_PLOT`, per subcohort and arm |
| `08_run_gistic.sh` | `RUN_GISTIC2`, per subcohort and arm |

The environment variables those scripts read (`SHEET_ID`, `TUMOUR_TYPE`, `THRESHOLD`, `GISTIC_INCLUDE_LIST`,
`GISTIC_EXCLUDE_LIST`, `EXCLUDE_FILE`) are now parameters: `cohort_prefix`, `recall_thresholds`,
`subcohorts.<name>.sample_list`, `subcohorts.<name>.exclude_list` and `exclude_samples`. `SHEET_ID` is gone -
the pipeline matches whatever directory fur-cnvkit names after the study id in the manifest.

## Inputs

- A tab-delimited metadata manifest describing the cohort. Three columns are read: the sample name
  (`Sanger_DNA_ID`, e.g. PD\*\*\*\*\*), the tumour/normal status (`Phenotype`, T/N) and the patient sex
  (`Sex`, M/F/U). Set `col_sample_id`, `col_TN` and `col_sex` if your manifest names them differently.

  ```{note}
  An `.xlsx` manifest is rejected. The pipeline reads the sex column itself, to tell `cnvkit call` whether to
  expect one or two sex chromosomes per sample. Export the sheet to TSV and point at the export.
  ```

- Indexed BAMs for every sample to be analysed, matched by `bam_path`.
- The list of tumours the study has agreed to analyse (`all_samples`), and one or more `subcohorts` to report
  on.
- A reference genome FASTA and its index, the baitset regions BED and a refFlat file.
- A pooled normal reference in CNVkit `.cnn` format, and the four static files it was built from.

Sensible Dermatlas defaults for everything except the cohort's own files are set in the `farm22` profile and
do not need changing unless you are using a non-standard reference genome or exome capture kit.

```{note}
The Dermatlas pooled reference was generated from 25 male and 25 female normals sequenced as part of the
hidradenoma cohort (6513_2706), on the WES5 baitset. A cohort on a different baitset or sequencing modality
needs its own reference - see [Building a pooled reference](#building-a-pooled-reference).
```

## Running the pipeline

See the repository README for the launcher, its environment contract and the reporting toggles. In the
managed case there is nothing to do beyond triggering the run from the
[Dermatlas cohorts page](https://team113.sanger.ac.uk/dermatlas/cohorts/).

## What each step does

### Parameter file

`CNVKIT_PARAMETER_FILE` runs `fur_cnvkit cnvkit_static_files --parameter-file-only`, which categorises the
cohort's BAMs by tumour/normal status and records them, along with the static files, in a `parameters.json`
that every later fur-cnvkit step reads.

The paths in that file are rewritten to bare filenames (`bin/rebase_parameter_paths.py`). fur-cnvkit resolves
them relative to the working directory, and under Nextflow each step runs in a directory of its own, so a
recorded absolute path would point at a task directory that no longer exists.

### Calling

`CNVKIT_CALL_COPY_NUMBER` runs the fur-cnvkit calling pipeline for the whole cohort at once, against the male
or female pooled reference according to each sample's sex. Per sample it produces coverage, ratio and segment
files, a scatter plot and a diagram; the file everything downstream uses is

```
<sample>.sample.dupmarked.call.thresholds_<thresholds>.call.median_centred.weight_ge_<weight>.0.cns
```

Two parameters shape this: `call_thresholds` (passed to `cnvkit call -t`) and `weight_filter_threshold`,
which drops segments whose `weight` is below 30 - the ones built from few or noisy bins, and the usual source
of spurious focal calls. `--skip-genemetrics` is always passed: Dermatlas does not use those outputs.

The scatter plots are worth looking at to judge how well a sample called at all.

### Hypersegmented samples

Some CNVkit profiles come back over-segmented - split into many small regions where the chromosomal mean
looks flat. They distort every GISTIC2 peak in a cohort, so `FIND_HYPERSEGMENTED_SAMPLES` reads the cohort's
own `*.metrics.tsv` and excludes any sample whose **average** segment is shorter than
`min_average_segment_size` (4 Mb). Two lists are published:

```
<prefix>_filtered_samples.txt            # kept
<prefix>_filtered_samples.excluded.txt   # hypersegmented, excluded from every subcohort
```

### Preliminary penetrance plots

`PRELIMINARY_PENETRANCE_PLOT` draws one genome-wide penetrance plot per purity in `purity_vs_log2ratio.tsv`,
over the retained samples. These are the QC plots the threshold choice is made from.

The thresholds are theoretical cut-offs:

```
log2(FC) = log2( (p * CN_tumour + (1 - p) * 2) / 2 )

# p = tumour purity as a fraction, CN_tumour = copy number in the tumour,
# and the denominator assumes a normal diploid copy number.
```

`purity_vs_log2ratio.tsv` tabulates that for CN=0 to 5 (with -infinity written as -1000000):

```
#Purity	CN0	CN1	CN2	CN3	CN4	CN5
0.30	-0.515	-0.234	0.000	0.202	0.379	0.536
0.40	-0.737	-0.322	0.000	0.263	0.485	0.678
0.50	-1.000	-0.415	0.000	0.322	0.585	0.807
...
```

Plots are named `penetrance_plot.purity<N>_gain<X>_loss<Y>.pdf`, with the binned counts beside them as
`.csv`. Gains are CN3 or higher and losses CN1 or lower, assuming tumour ploidy 2.

What to look for:

- many small recurrent gains/losses near centromeres, telomeres, reference gaps or other difficult regions
  are more likely artefacts than the large calls
- if raising the thresholds costs you chromosome-arm-level events, they are already too stringent
- these are for QC; further filtering happens per subcohort below

When a different purity fits the cohort better, set `recall_thresholds` to that row's CN0, CN1, CN3 and CN4
values and run again.

### Re-calling

`CNVKIT_RECALL` runs `cnvkit call --method none -m threshold -t=<recall_thresholds>` per sample, adding `-y`
for male samples so their single-copy X and Y are not called as losses. `--method none` leaves the log2 ratios
untouched: only the integer `cn` column changes. When `recall_thresholds` equals `call_thresholds` the values
are unchanged, which is why the default is a no-op rather than a special case.

Samples whose manifest sex is neither M nor F are called but not re-called, and a warning names them.

### Filtering and GISTIC2 inputs

For each subcohort, `FILTER_SEGMENTS` runs twice.

The **filtered** arm drops segments that overlap `difficult_regions` by at least `bed_overlap_threshold`
(60%) of their length, that are shorter than `min_segment_size` (400 kb), or whose weight is below
`segment_weight`. It publishes the surviving segments, a `*.removed.cns` beside each one listing what was
taken out, and the GISTIC2 `.seg` export.

The **unfiltered** arm exports the same samples with nothing dropped.

```{note}
Read the `*.removed.cns` files. If a known event - a focal *MYC* amplification, a *TP53* deletion - has been
filtered out, the thresholds are too aggressive for this cohort.
```

Both arms then get a penetrance plot, drawn over exactly the samples in their `.seg` file, so the effect of
filtering is visible side by side.

### GISTIC2

`RUN_GISTIC2` scores each arm of each subcohort, with `-ta`/`-td` taken from the same thresholds used to
re-call `cn` (as positive amplitudes). Results are published under
`<subcohort>/gistic2/{filtered,unfiltered}/`.

The filtered results should generally be the less noisy of the two; the unfiltered ones are there to show
whether over-filtering has occurred.

### GISTIC2 QC

Further QC of the GISTIC2 result - comparing focal peaks back to the CNVkit segments, assessing overlap with
difficult regions, and annotating genes and cancer genes within peaks - is done with the Dermatlas
[`gistic_assess`](https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/gistic_assess) scripts, and
is not part of this pipeline. Run them over both arms:

```bash
IFS=',' read -ra v <<< "${THRESHOLD}"
GISTIC_QC_THRESH="${v[1]},${v[2]}"

for type in filtered unfiltered; do
  bsub -o logs/gistic_assess.${type}.%J.o -e logs/gistic_assess.${type}.%J.e \
    -M4000 -R"select[mem>4000] rusage[mem=4000]" \
    "Rscript ${PROJECT_DIR}/scripts/gistic_assess/gistic2_filter.R \
      --prefix ${PREFIX} \
      --gistic-all-lesions-file <outdir>/<subcohort>/gistic2/${type}/all_lesions.conf_95.txt \
      --gistic-segments-file <outdir>/<subcohort>/${type}/<prefix>_${type}_segments.seg \
      --residual-q-value-cutoff 0.1 \
      --output-dir <outdir>/<subcohort>/gistic2/${type}/ \
      -d ${RESOURCES_DIR}/gistic/gistic2_difficult_regions.bed \
      -e ${RESOURCES_DIR}/ensembl/Homo_sapiens.GRCh38.103.chr.gff3.gz \
      -c ${RESOURCES_DIR}/gistic/cancer_gene_matrix.tsv \
      --oncokb ${RESOURCES_DIR}/oncokb/cancerGeneList.nocosmic.tsv \
      --gistic-log2=${GISTIC_QC_THRESH}"
done
```

```{note}
Check the logs for the parameters actually used. The sample/CN concordance threshold between GISTIC2 and
CNVkit is 0.5, against 0.75 for ASCAT: GISTIC2 and CNVkit need not agree closely here, because a single log2
ratio threshold is applied across the whole cohort (there being no good per-sample purity/ploidy estimate),
whereas ASCAT fits purity and ploidy per sample.
```

To see the calls that pass QC:

```bash
awk -F "\t" '$13~/PASS|qc_status/' <subcohort>/gistic2/filtered/<prefix>_gistic_cohort_summary.tsv
```

## Building a pooled reference

`-entry BUILD_REFERENCE` generates the static files and the pooled normal reference for a baitset that does
not have one. It is a long-running, all-vs-all comparison of the cohort's normals - roughly four days for 60
WES samples on 60 CPUs - and is run once per baitset, not once per cohort.

```bash
nextflow run main.nf -entry BUILD_REFERENCE -profile farm22 \
  --bam_path "/path/to/bams/**bam{,.bai}" \
  --metadata_manifest /path/to/metadata.tsv \
  --baitset_bed /path/to/new_baitset.bed \
  --outdir /path/to/reference_output
```

It publishes:

```
static_files/
├── access-<genome>.bed
├── <baitset>.target.bed
├── <baitset>.antitarget.bed
├── <baitset>.baitset_genes.txt
└── parameters.json
pooled_reference/
├── female/{female.reference.cnn,samples_used_in_reference.txt,normal_vs_normal/}
└── male/{male.reference.cnn,samples_used_in_reference.txt,normal_vs_normal/}
```

The `.cnn` files are built from the least noisy 80% of the normals of each sex, as judged by normal-vs-normal
comparison; `samples_used_in_reference.txt` records which ones. Point `access_bed`, `targets_bed`,
`antitargets_bed`, `baitset_genes_file`, `male_reference` and `female_reference` at these outputs for every
cohort on that baitset - the reference is indexed by the bins in those BEDs, so the two sets must always
travel together.

`sv_blacklist` is subtracted from the accessible regions. Dermatlas uses a merged file of the 10x Genomics
`sv_blacklist.bed` for GRCh38 and the UCSC "Unusual Regions" track from the Problematic Regions set.
