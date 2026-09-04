include { CNVKIT_STATIC_FILES; CNVKIT_POOLED_REFERENCE } from '../../modules/cnvkit.nf'

// Build the static files and the pooled normal reference a cohort is called against.
//
// This is the alternative entry point (`-entry BUILD_REFERENCE`) and is deliberately not
// part of a calling run: an all-vs-all comparison of the cohort's normals takes days,
// and the resulting .cnn files are then reused by every cohort on the same baitset. Run
// it once per baitset, publish the two .cnn files, and point `male_reference` /
// `female_reference` at them for every run thereafter.
workflow BUILD_REFERENCE {
    take:
    cohort_bams          // value channel of tuple(meta, [bams], [bais])
    metadata_manifest
    reference_genome
    reference_genome_index
    baitset_bed
    refflat_file
    sv_blacklist
    exclude_samples

    main:

    CNVKIT_STATIC_FILES(
        cohort_bams,
        metadata_manifest,
        reference_genome,
        reference_genome_index,
        baitset_bed,
        refflat_file,
        sv_blacklist,
        exclude_samples
    )

    // The BED files and the baitset gene list are named in the parameter file by
    // basename, so they have to be staged alongside it. They are collected into a single
    // `static_files` input rather than four named ones because the reference builder
    // reads them through the parameter file, never by argument.
    //
    // CNVKIT_STATIC_FILES takes only value channels, so it runs once and its outputs are
    // value channels too - which is what lets `parameters` be both consumed here and
    // emitted below. The joined `static_files` is a queue channel, hence the .first().
    parameters = CNVKIT_STATIC_FILES.out.parameters
    static_files = CNVKIT_STATIC_FILES.out.access
        .join(CNVKIT_STATIC_FILES.out.targets)
        .join(CNVKIT_STATIC_FILES.out.antitargets)
        .join(CNVKIT_STATIC_FILES.out.baitset_genes)
        .map { meta, access, targets, antitargets, baitset_genes ->
            [access, targets, antitargets, baitset_genes]
        }
        .first()

    reference_input = parameters
        .combine(cohort_bams)
        .map { meta, parameter_file, _bam_meta, bams, bais -> tuple(meta, parameter_file, bams, bais) }

    CNVKIT_POOLED_REFERENCE(
        reference_input,
        metadata_manifest,
        reference_genome,
        reference_genome_index,
        baitset_bed,
        static_files
    )

    emit:
    parameters       = parameters
    // The four static files as one list, in the order access, targets, antitargets,
    // baitset genes. They are published by CNVKIT_STATIC_FILES regardless.
    static_files     = static_files
    male_reference   = CNVKIT_POOLED_REFERENCE.out.male_reference
    female_reference = CNVKIT_POOLED_REFERENCE.out.female_reference
}
