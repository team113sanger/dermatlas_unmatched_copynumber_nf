#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { INGEST_COHORT } from './subworkflows/local/ingest_cohort.nf'
include { CALL_COPY_NUMBER } from './subworkflows/local/call_copy_number.nf'
include { ANALYSE_SUBCOHORT } from './subworkflows/local/analyse_subcohort.nf'
include { BUILD_REFERENCE as BUILD_POOLED_REFERENCE } from './subworkflows/local/build_reference.nf'

include { countDataRows; readHeader; parseThresholds } from './modules/utils.nf'

// Sentinel for the optional file inputs. Nextflow has no "no file" for a `path` input, so
// an empty placeholder is staged and the processes test for its name.
def noFile() {
    return file("${projectDir}/assets/NO_FILE", checkIfExists: true)
}

// Validate everything both entry points need, and return the values they share.
//
// All of this runs before a single task is scheduled. That matters more here than in most
// pipelines: the calling step is hours of compute and the reference builder is days of it,
// so a mistyped column name or an empty sample list has to be caught now, not by an
// output glob that matches nothing at the end.
def validateCohortInputs() {
    if (!params.bam_path) {
        error("params.bam_path must be set: a glob matching the cohort's indexed BAMs and their " +
              "indexes, e.g. \"/path/to/bams/**bam{,.bai}\".")
    }
    if (!params.metadata_manifest) {
        error("params.metadata_manifest must be set: the tab-delimited metadata manifest for the cohort.")
    }

    def manifest = file(params.metadata_manifest, checkIfExists: true)

    // fur-cnvkit reads .xlsx manifests, but this pipeline also reads the sex column
    // itself - CNVKIT_RECALL has to be told per sample whether to expect one or two sex
    // chromosomes - and a spreadsheet cannot be parsed from the workflow.
    if (manifest.name.toLowerCase().endsWith('.xlsx')) {
        error("params.metadata_manifest '${params.metadata_manifest}' is an .xlsx file. This pipeline " +
              "needs a tab-delimited manifest; export the sheet to TSV and point at that instead.")
    }
    if (countDataRows(manifest, 1) == 0) {
        error("No sample records found in params.metadata_manifest '${params.metadata_manifest}'. " +
              "Expected a header row followed by one row per sample.")
    }

    // A mistyped column name is indistinguishable from an empty manifest once the rows
    // have been parsed, so the header is checked by name.
    def header = readHeader(manifest, "\t")
    def required = [params.col_sample_id, params.col_sex, params.col_TN]
    def missing = required.findAll { column -> !header.contains(column) }
    if (missing) {
        error("Metadata manifest '${params.metadata_manifest}' is missing required column(s): " +
              "${missing.join(', ')}. Columns found: ${header.join(', ')}. Set --col_sample_id / " +
              "--col_sex / --col_TN to match your manifest.")
    }

    // Reference data. These have site defaults in the farm22 profile, so an unset one
    // means the run is on another profile and has not been told where its references are.
    def reference_params = ['reference_genome', 'reference_genome_index', 'baitset_bed', 'refflat_file']
    def unset_references = reference_params.findAll { name -> !params[name] }
    if (unset_references) {
        error("params.${unset_references.join(', params.')} must be set. The farm22 profile supplies " +
              "the Dermatlas defaults; on any other profile they have to be given explicitly.")
    }

    // The ids named here are dropped before anything is scheduled, so an excluded sample
    // costs no compute at all.
    def excluded_ids = []
    if (params.exclude_samples) {
        excluded_ids = file(params.exclude_samples, checkIfExists: true)
            .readLines()
            .collect { line -> line.trim().split('\t')[0].trim() }
            .findAll { id -> id }
        log.info("Excluding ${excluded_ids.size()} sample(s) listed in ${params.exclude_samples}")
    }

    return [
        manifest: manifest,
        excluded_ids: excluded_ids,
        exclude_file: params.exclude_samples ? file(params.exclude_samples, checkIfExists: true) : noFile(),
        reference_genome: file(params.reference_genome, checkIfExists: true),
        reference_genome_index: file(params.reference_genome_index, checkIfExists: true),
        baitset_bed: file(params.baitset_bed, checkIfExists: true),
        refflat_file: file(params.refflat_file, checkIfExists: true)
    ]
}

// Default entry point: call a cohort against an existing pooled normal reference.
workflow UNMATCHED_COPY_NUMBER {

    def inputs = validateCohortInputs()

    if (!params.cohort_prefix) {
        error("params.cohort_prefix must be set. It prefixes every sample list, .seg file and plot " +
              "this run produces, so leaving it unset would publish files named 'null_...'.")
    }
    if (!params.all_samples) {
        error("params.all_samples must be set: the list of tumour samples the study has agreed to " +
              "analyse. It decides which BAMs are called at all.")
    }
    if (!params.male_reference || !params.female_reference) {
        error("params.male_reference and params.female_reference must both be set. To build them from " +
              "this cohort's own normals instead, run with `-entry BUILD_REFERENCE`.")
    }

    // The four static files must be the ones the pooled reference was built from: the
    // .cnn reference is indexed by the bins in these BEDs, so calling against a reference
    // built from a different target set silently compares incomparable bins.
    def static_file_params = ['access_bed', 'targets_bed', 'antitargets_bed', 'baitset_genes_file',
                              'purity_log2_table', 'difficult_regions', 'gistic_refgene_file']
    def unset_static_files = static_file_params.findAll { name -> !params[name] }
    if (unset_static_files) {
        error("params.${unset_static_files.join(', params.')} must be set. The four static files must " +
              "be the ones generated alongside the pooled reference named by params.male_reference / " +
              "params.female_reference; the farm22 profile supplies the Dermatlas defaults for all of these.")
    }

    // The thresholds the run re-calls `cn` with. They default to the ones CNVkit called
    // with, in which case re-calling is a no-op that leaves the calls unchanged - the
    // preliminary penetrance plots are what tell you whether a different purity fits the
    // cohort better.
    def thresholds = parseThresholds(params.recall_thresholds ?: params.call_thresholds, 'params.recall_thresholds')
    log.info("Calling thresholds: ${params.call_thresholds}; re-calling with: ${thresholds.thresholds} " +
             "(loss ${thresholds.loss}, gain ${thresholds.gain})")

    if (!params.subcohorts) {
        error("params.subcohorts must define at least one subcohort, e.g.\n" +
              "  subcohorts = ['one_per_patient': [sample_list: '/path/to/samples.tsv']]")
    }

    // Reject malformed subcohorts by name, before they become a null path in a channel.
    def malformed = params.subcohorts.findAll { _name, config ->
        !(config instanceof Map) || !config.sample_list
    }
    if (malformed) {
        error("Malformed 'subcohorts' entries: ${malformed.keySet().join(', ')}. Each subcohort must " +
              "define a 'sample_list'.")
    }

    // An empty sample list would produce an empty .seg file and a GISTIC2 run with nothing
    // to score, so a subcohort that has lost all of its samples is dropped with a warning
    // rather than failing the run for the subcohorts that still have some.
    def subcohort_sets = params.subcohorts.collect { name, config ->
        def sample_list = file(config.sample_list, checkIfExists: true)
        if (countDataRows(sample_list, 0) == 0) {
            log.warn("Skipping subcohort '${name}': sample list '${config.sample_list}' contains no samples.")
            return null
        }
        // An empty string, not null: a null inside a channel item is not carried reliably.
        return tuple(name, sample_list, config.exclude_list ?: '')
    }.findAll { entry -> entry != null }

    if (!subcohort_sets) {
        error("All configured subcohorts (${params.subcohorts.keySet().join(', ')}) have empty sample lists.")
    }
    log.info("Analysing subcohorts: ${subcohort_sets.collect { name, _list, _exclude -> name }.join(', ')}")

    cohort_meta = ["study_id": params.study_id, "prefix": params.cohort_prefix]

    INGEST_COHORT(
        params.bam_path,
        inputs.manifest,
        params.all_samples,
        inputs.excluded_ids,
        cohort_meta
    )

    CALL_COPY_NUMBER(
        INGEST_COHORT.out.cohort_bams,
        INGEST_COHORT.out.sample_sex,
        inputs.manifest,
        inputs.reference_genome,
        inputs.reference_genome_index,
        inputs.baitset_bed,
        inputs.refflat_file,
        file(params.access_bed, checkIfExists: true),
        file(params.targets_bed, checkIfExists: true),
        file(params.antitargets_bed, checkIfExists: true),
        file(params.baitset_genes_file, checkIfExists: true),
        file(params.male_reference, checkIfExists: true),
        file(params.female_reference, checkIfExists: true),
        file(params.all_samples, checkIfExists: true),
        inputs.exclude_file,
        params.purity_log2_table,
        thresholds
    )

    // A subcohort with no exclude list of its own excludes the samples this run found to
    // be hypersegmented. That is the default deliberately: those samples are noise in any
    // subcohort they appear in, and the list is only known once calling has finished.
    subcohorts_ch = Channel.fromList(subcohort_sets)
        .combine(CALL_COPY_NUMBER.out.hypersegmented)
        .map { name, sample_list, exclude_list, _meta, hypersegmented ->
            tuple(
                cohort_meta + ["cohort_id": name],
                sample_list,
                exclude_list ? file(exclude_list, checkIfExists: true) : hypersegmented
            )
        }

    ANALYSE_SUBCOHORT(
        subcohorts_ch,
        CALL_COPY_NUMBER.out.recalled_cns,
        file(params.difficult_regions, checkIfExists: true),
        file(params.gistic_refgene_file, checkIfExists: true),
        thresholds
    )
}

// Alternative entry point: build the static files and pooled normal reference this
// baitset is called against. Run once per baitset, not once per cohort.
workflow BUILD_REFERENCE {

    def inputs = validateCohortInputs()

    if (!params.sv_blacklist) {
        error("params.sv_blacklist must be set when building a reference: the BED of regions to keep " +
              "out of the accessible-regions file (blacklisted and problematic regions).")
    }

    INGEST_COHORT(
        params.bam_path,
        inputs.manifest,
        // No sample universe: the reference is built from the cohort's normals, which are
        // by definition absent from the list of tumours to analyse.
        null,
        inputs.excluded_ids,
        ["study_id": params.study_id, "prefix": params.cohort_prefix ?: params.study_id]
    )

    BUILD_POOLED_REFERENCE(
        INGEST_COHORT.out.cohort_bams,
        inputs.manifest,
        inputs.reference_genome,
        inputs.reference_genome_index,
        inputs.baitset_bed,
        inputs.refflat_file,
        file(params.sv_blacklist, checkIfExists: true),
        inputs.exclude_file
    )
}

workflow {
    UNMATCHED_COPY_NUMBER()
}

workflow.onComplete {
    // Runs on both success and failure, after all processes have finished.
    // All reporting (Slack + analysis-log) is handled in one reusable call.
    Utils.reportRun(workflow, params)
}
