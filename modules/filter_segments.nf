// Segment-level filtering, driven by bin/filter_cnvkit.py.
//
// The script has two distinct jobs, split here into two processes because they run at
// different points and over different inputs:
//   - FIND_HYPERSEGMENTED_SAMPLES reads the cohort's CNVkit metrics and decides which
//     samples are too finely segmented to be usable at all.
//   - FILTER_SEGMENTS takes the re-called .cns files for one subcohort and produces the
//     GISTIC2 input .seg file, optionally dropping segments that overlap difficult
//     regions or fall below the size/weight cut-offs.

process FIND_HYPERSEGMENTED_SAMPLES {
    label 'process_tiny'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    publishDir path: "${params.outdir}/cnvkit_cn_calling",
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(cns_files), path(metrics)
    path(all_samples)

    output:
    tuple val(meta), path("${meta.prefix}_filtered_samples.txt"), emit: retained
    tuple val(meta), path("${meta.prefix}_filtered_samples.excluded.txt"), emit: excluded

    script:
    // --min-seg-size is a per-sample *average* segment size: a sample whose mean segment
    // is shorter than this is treated as hypersegmented and dropped. It is unrelated to
    // FILTER_SEGMENTS' --remove-seg, which drops individual short segments.
    """
    filter_cnvkit.py \\
        --input ${cns_files} \\
        --samples ${all_samples} \\
        --quality ${metrics} \\
        --min-seg-size ${params.min_average_segment_size} \\
        --sample-list-out ${meta.prefix}_filtered_samples.txt
    """

    stub:
    """
    cut -f1 ${all_samples} > ${meta.prefix}_filtered_samples.txt
    touch ${meta.prefix}_filtered_samples.excluded.txt
    """
}

process FILTER_SEGMENTS {
    label 'process_low'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    publishDir path: { "${params.outdir}/${meta.cohort_id}/${meta.analysis_type}" },
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(cns_files), path(include_samples), path(exclude_samples)
    path(difficult_regions)

    output:
    // Only the 'filtered' arm rewrites the .cns files; 'unfiltered' exports the inputs
    // to GISTIC2 untouched, so it produces no per-sample .cns of its own.
    tuple val(meta), path("cns/*.${meta.analysis_type}.cns"), emit: cns, optional: true
    tuple val(meta), path("cns/*.removed.cns"), emit: removed, optional: true
    tuple val(meta), path("${meta.prefix}_${meta.analysis_type}_segments.seg"), emit: seg
    tuple val(meta), path("${meta.prefix}_gistic_${meta.analysis_type}_samples.txt"), emit: samples

    script:
    // 'unfiltered' keeps every segment and exists purely as the comparison arm: the same
    // samples, the same .seg export, but no segment dropped. Running GISTIC2 over both is
    // how over-filtering is detected.
    //
    // The argument list is assembled in Groovy rather than interpolated line by line:
    // an empty interpolation on its own line would end the shell command at the
    // preceding backslash and silently drop every argument after it.
    def filter_arguments = []
    filter_arguments << "--input ${cns_files}"
    filter_arguments << "--samples ${include_samples}"
    if (exclude_samples.name != 'NO_FILE') {
        filter_arguments << "--exclude-samples ${exclude_samples}"
    }
    if (meta.analysis_type == 'filtered') {
        filter_arguments << "--filter-cns"
        filter_arguments << "--weight ${params.segment_weight}"
        filter_arguments << "--bed-file ${difficult_regions}"
        filter_arguments << "--bed-threshold ${params.bed_overlap_threshold}"
        filter_arguments << "--remove-seg ${params.min_segment_size}"
    }
    filter_arguments << "--suffix ${meta.analysis_type}"
    filter_arguments << "--output-dir cns"
    filter_arguments << "--gistic-seg ${meta.prefix}_${meta.analysis_type}_segments.seg"
    filter_arguments << "--strip-filenames"
    filter_arguments << "--sample-list-out ${meta.prefix}_gistic_${meta.analysis_type}_samples.txt"
    """
    mkdir -p cns
    filter_cnvkit.py \\
        ${filter_arguments.join(' \\\n        ')}

    # Shorten the names CNVkit accumulated one step at a time, e.g.
    #   PD60062a.sample.dupmarked.call.thresholds_X.call.median_centred.weight_ge_30.0.filtered.cns
    #   -> PD60062a.thresholds_X.med_centred.wt_ge_30.0.filtered.cns
    for file in cns/*.cns; do
        new_name="\${file//sample.dupmarked.call./}"
        new_name="\${new_name//call.median_centred.weight/med_centred.wt}"
        [ "\${file}" = "\${new_name}" ] || mv "\${file}" "\${new_name}"
    done
    """

    stub:
    def write_cns = meta.analysis_type == 'filtered'
    // Honour the exclude list, as the script does, so a subcohort that loses every sample
    // can be exercised without the real tools.
    def drop_excluded = exclude_samples.name != 'NO_FILE' ? "| { grep -vxFf <(cut -f1 ${exclude_samples}) || true; }" : ''
    """
    mkdir -p cns
    if ${write_cns}; then
        for cns in ${cns_files}; do
            sample=\$(echo \${cns} | cut -d. -f1)
            echo stub > cns/\${sample}.${meta.analysis_type}.cns
            echo stub > cns/\${sample}.removed.cns
        done
    fi
    echo stub > ${meta.prefix}_${meta.analysis_type}_segments.seg
    cut -f1 ${include_samples} ${drop_excluded} > ${meta.prefix}_gistic_${meta.analysis_type}_samples.txt
    """
}
