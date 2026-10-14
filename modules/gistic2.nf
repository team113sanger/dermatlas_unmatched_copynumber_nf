// GISTIC2 over the per-subcohort segment files.
//
// One instance runs per (subcohort, analysis_type) pair, so the filtered and unfiltered
// segments are scored independently and their peak calls can be compared. GISTIC2 writes
// a fixed set of filenames into `-b`, which is the task directory here.

process RUN_GISTIC2 {
    label 'process_medium'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/analysis-methods/gistic2:0.5.0"
    publishDir path: { "${params.outdir}/${meta.cohort_id}/gistic2/${meta.analysis_type}" },
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(segment_file)
    path(refgene_file)

    // Outputs are optional so a subcohort GISTIC2 cannot analyse (see script) is
    // skipped downstream rather than failing the pipeline
    output:
    tuple val(meta), path("all_lesions.conf_95.txt"), emit: lesions, optional: true
    tuple val(meta), path("broad_significance_results.txt"), emit: broad, optional: true
    tuple val(meta), path("broad_values_by_arm.txt"), emit: arms, optional: true
    tuple val(meta), path("*.gistic"), emit: scores, optional: true
    tuple val(meta), path("*.txt"), emit: tables, optional: true
    tuple val(meta), path("*.pdf"), emit: pdfs, optional: true
    tuple val(meta), path("*.png"), emit: plots, optional: true
    tuple val(meta), path("*.mat"), emit: mats, optional: true
    tuple val(meta), path("gistic2_skipped.log"), emit: skipped, optional: true

    script:
    // -ta / -td are amplitude thresholds and both are given as positive numbers, so the
    // loss threshold - a negative log2 ratio everywhere else in this pipeline - is
    // negated here. They come from the same threshold string used to re-call `cn`, so
    // GISTIC2 scores exactly the gains and losses that CNVkit called.
    """
    set +e
    /opt/repo/gp_gistic2_from_seg \\
        -b . \\
        -seg ${segment_file} \\
        -refgene ${refgene_file} \\
        -genegistic 1 \\
        -smallmem 1 \\
        -broad 1 \\
        -brlen ${params.gistic_broad_length_cutoff} \\
        -conf ${params.gistic_confidence_level} \\
        -armpeel 1 \\
        -savegene 1 \\
        -gcm extreme \\
        -v 20 \\
        -ta ${meta.gain} \\
        -td ${meta.loss_amplitude} 2>&1 | tee gistic2.log
    gistic_status=\${PIPESTATUS[0]}
    set -e

    # GISTIC2 needs >1 sample with copy number signal; a subcohort where every
    # sample but one is flat is not an error, just nothing to analyse
    if grep -q "Only one sample remains after noise filtering" gistic2.log; then
        msg="WARNING: GISTIC2 skipped for ${meta.cohort_id} (${meta.analysis_type}): only one sample remains after noise filtering."
        echo "\$msg" | tee gistic2_skipped.log >&2
        exit 0
    fi

    exit \$gistic_status
    """

    stub:
    """
    echo stub > all_lesions.conf_95.txt
    echo stub > broad_significance_results.txt
    echo stub > broad_values_by_arm.txt
    echo stub > scores.gistic
    """
}
