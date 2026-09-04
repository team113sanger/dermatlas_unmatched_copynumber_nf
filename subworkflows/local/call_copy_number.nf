include { CNVKIT_PARAMETER_FILE; CNVKIT_CALL_COPY_NUMBER; CNVKIT_RECALL } from '../../modules/cnvkit.nf'
include { FIND_HYPERSEGMENTED_SAMPLES } from '../../modules/filter_segments.nf'
include { PENETRANCE_PLOT as PRELIMINARY_PENETRANCE_PLOT } from '../../modules/penetrance.nf'

// Call copy number for a whole cohort against a pre-computed pooled normal reference,
// then decide which samples are usable and what the `cn` calls should be.
//
// This is the part of the original workflow that a human drove step by step: call, look
// at a sweep of penetrance plots, pick a purity, re-call with its thresholds. The sweep
// is still produced - one plot per purity, in parallel - but the thresholds the run
// re-calls with come from params.recall_thresholds, so the whole thing completes
// unattended and the plots are there to justify (or revise) that choice next time.
//
// fur-cnvkit takes the whole cohort in one invocation, so most of the channels below
// carry a single cohort-level item. Those are turned into value channels with .first()
// as soon as they come out of a process: a queue channel can only be consumed once, and
// each of these feeds two or three consumers.
workflow CALL_COPY_NUMBER {
    take:
    cohort_bams          // value channel of tuple(meta, [bams], [bais])
    sample_sex           // tuple(sample_id, sex) per tumour
    metadata_manifest
    reference_genome
    reference_genome_index
    baitset_bed
    refflat_file
    access_bed
    targets_bed
    antitargets_bed
    baitset_genes_file
    male_reference
    female_reference
    all_samples
    exclude_samples
    purity_log2_table
    recall_thresholds    // parsed threshold map (see modules/utils.nf)

    main:

    CNVKIT_PARAMETER_FILE(
        cohort_bams,
        metadata_manifest,
        reference_genome,
        reference_genome_index,
        baitset_bed,
        refflat_file,
        access_bed,
        targets_bed,
        antitargets_bed,
        baitset_genes_file,
        exclude_samples
    )

    // The parameter file records its BAMs by basename, so the calling process has to
    // stage exactly the same set of files for those names to resolve.
    calling_input = CNVKIT_PARAMETER_FILE.out.parameters
        .combine(cohort_bams)
        .map { meta, parameters, _bam_meta, bams, bais -> tuple(meta, parameters, bams, bais) }

    CNVKIT_CALL_COPY_NUMBER(
        calling_input,
        metadata_manifest,
        baitset_genes_file,
        male_reference,
        female_reference
    )

    called_cns = CNVKIT_CALL_COPY_NUMBER.out.cns.first()
    called_metrics = CNVKIT_CALL_COPY_NUMBER.out.metrics.first()

    // Hypersegmented samples are identified from the cohort's own metrics rather than
    // from a fixed list: a sample whose average segment is short has a profile too noisy
    // for GISTIC2, and including it distorts every peak in the cohort.
    FIND_HYPERSEGMENTED_SAMPLES(
        called_cns
            .combine(called_metrics)
            .map { meta, cns, _metrics_meta, metrics -> tuple(meta, cns, metrics) },
        all_samples
    )

    retained_samples = FIND_HYPERSEGMENTED_SAMPLES.out.retained.first()
    hypersegmented_samples = FIND_HYPERSEGMENTED_SAMPLES.out.excluded.first()

    // The theoretical log2 ratio thresholds for CN=1 and CN=3 at each tumour purity. The
    // table's CN2 column is always 0.0 (a diploid genome is the reference), so only the
    // CN1 (loss) and CN3 (gain) columns are read.
    purity_sweep = Channel.fromPath(purity_log2_table, checkIfExists: true)
        .splitCsv(sep: '\t')
        .filter { row -> row && !row[0].startsWith('#') }
        .map { row ->
            [
                "label": "purity${row[0]}",
                "gain": row[4],
                "loss": row[2],
                "plot_dir": "cnvkit_cn_penetrance_plots"
            ]
        }

    // Every plot in the sweep is drawn from the same segments and the same sample list;
    // only the gain/loss thresholds differ, so the sweep is a plain combine().
    preliminary_plot_input = purity_sweep
        .combine(called_cns)
        .combine(retained_samples)
        .map { plot_meta, _cns_meta, cns, _retained_meta, retained ->
            tuple(plot_meta, cns, retained, exclude_samples)
        }

    PRELIMINARY_PENETRANCE_PLOT(preliminary_plot_input)

    // Re-call `cn` per sample with the chosen thresholds. `cnvkit call --method none`
    // leaves the log2 ratios alone, so this only rewrites the integer copy-number column
    // - but it has to be told each sample's sex, which is why it is per-sample here while
    // the calling step above was per-cohort.
    recall_input = called_cns
        .flatMap { meta, cns ->
            [cns].flatten().collect { segments -> tuple(segments.name.tokenize('.').first(), segments) }
        }
        .join(sample_sex)
        .map { sample_id, cns, sex -> tuple(["sample_id": sample_id, "sex": sex], cns) }

    CNVKIT_RECALL(recall_input, recall_thresholds.thresholds)

    emit:
    parameters     = CNVKIT_PARAMETER_FILE.out.parameters
    cns            = called_cns
    metrics        = called_metrics
    retained       = retained_samples
    hypersegmented = hypersegmented_samples
    // Back to one bundle per cohort for the subcohort analyses.
    recalled_cns   = CNVKIT_RECALL.out.cns.map { meta, cns -> cns }.collect()
}
