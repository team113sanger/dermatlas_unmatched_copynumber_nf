// Genome-wide copy-number penetrance plots, drawn by bin/visualise_penetrance.py.
//
// The same process serves both places the plots are needed:
//   - the preliminary sweep, one instance per purity in the purity/log2-ratio table,
//     used to choose the gain/loss thresholds for the run;
//   - the final per-subcohort plots, drawn once each from the filtered and the
//     unfiltered segments so over-filtering is visible side by side.
//
// The script takes a directory and a glob rather than a file list, so the .cns files are
// staged into the task directory and read back out of it with `-d .`.

process PENETRANCE_PLOT {
    label 'process_low'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    publishDir path: { "${params.outdir}/${meta.plot_dir}" },
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(cns_files, stageAs: 'segments/*'), path(include_samples), path(exclude_samples)

    output:
    tuple val(meta), path("*.pdf"), emit: plot
    tuple val(meta), path("*.csv"), emit: plot_data

    script:
    // meta.label distinguishes the arms that share a plot directory (e.g. "purity0.40"
    // for the sweep, "filtered"/"unfiltered" for the final plots), so the two never
    // collide when publishDir copies them into the same place.
    def exclude_arg = exclude_samples.name != 'NO_FILE' ? "--exclude-samples ${exclude_samples}" : ''
    def name = "penetrance_plot.${meta.label}_gain${meta.gain}_loss${meta.loss}"
    """
    visualise_penetrance.py \\
        -d segments \\
        --pattern "*.cns" \\
        --samples ${include_samples} \\
        ${exclude_arg} \\
        --gain-threshold ${meta.gain} \\
        --loss-threshold ${meta.loss} \\
        --y-range ${params.penetrance_y_range} \\
        -o ${name}.pdf \\
        --save-data ${name}.csv \\
        --verbose
    """

    stub:
    def name = "penetrance_plot.${meta.label}_gain${meta.gain}_loss${meta.loss}"
    """
    echo stub > ${name}.pdf
    echo stub > ${name}.csv
    """
}
