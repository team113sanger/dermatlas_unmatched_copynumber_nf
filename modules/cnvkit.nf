// CNVkit steps, all run through the fur-cnvkit wrapper.
//
// Every process writes into its own task directory (`-o .`) and lets publishDir place
// the results, so the output layout is owned by this pipeline rather than by the
// absolute paths that the original wrapper scripts baked into each command.
//
// fur-cnvkit names its output directory after the study id it reads from the metadata
// manifest ("default" for a TSV manifest, the sheet name for an .xlsx one), so the
// output globs below match `*/` rather than a fixed directory: the pipeline never has
// to be told which of the two it is.

process CNVKIT_STATIC_FILES {
    label 'process_low'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    publishDir path: "${params.outdir}/static_files",
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(bams), path(bais)
    path(metadata_manifest)
    path(reference_genome)
    path(reference_genome_index)
    path(baitset_bed)
    path(refflat_file)
    path(sv_blacklist)
    path(exclude_samples)

    output:
    tuple val(meta), path("parameters.json"), emit: parameters
    tuple val(meta), path("*.target.bed"), emit: targets
    tuple val(meta), path("*.antitarget.bed"), emit: antitargets
    tuple val(meta), path("access-*.bed"), emit: access
    tuple val(meta), path("*.baitset_genes.txt"), emit: baitset_genes

    script:
    // `exclude_samples` is an optional input: Nextflow stages the empty list as no file
    // at all, so the flag is only added when a real file was provided.
    def exclude_arg = exclude_samples.name != 'NO_FILE' ? "-e ${exclude_samples}" : ''
    def contig_args = params.unplaced_contig_prefixes.collect { prefix -> "\"${prefix}\"" }.join(' ')
    """
    fur_cnvkit cnvkit_static_files \\
        -b ${bams} \\
        ${exclude_arg} \\
        -f ${reference_genome} \\
        -t ${baitset_bed} \\
        -u ${contig_args} \\
        -x ${sv_blacklist} \\
        -r ${refflat_file} \\
        -m ${metadata_manifest} \\
        -o . \\
        --metadata-sample-id-column "${params.col_sample_id}" \\
        --metadata-tumour-normal-column "${params.col_TN}" \\
        --metadata-sex-column "${params.col_sex}" \\
        --verbose DEBUG

    rebase_parameter_paths.py --input parameters.json --output parameters.rebased.json
    mv parameters.rebased.json parameters.json
    """

    stub:
    """
    echo stub > baitset.target.bed
    echo stub > baitset.antitarget.bed
    echo stub > access-genome.bed
    echo stub > baitset.baitset_genes.txt
    cat <<'JSON' > parameters.json
    {
        "all_bams": ["stub.bam"],
        "tumour_bams": ["stub.bam"],
        "normal_bams": [],
        "reference_fasta": "genome.fa",
        "unplaced_contig_prefixes": [],
        "baitset_bed": "baitset.bed",
        "refflat_file": "refFlat.txt",
        "sample_metadata_xlsx": "metadata.tsv",
        "access_bed": "access-genome.bed",
        "targets_bed": "baitset.target.bed",
        "antitargets_bed": "baitset.antitarget.bed",
        "baitset_genes_file": "baitset.baitset_genes.txt"
    }
    JSON
    """
}

process CNVKIT_PARAMETER_FILE {
    label 'process_low'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    publishDir path: "${params.outdir}",
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(bams), path(bais)
    path(metadata_manifest)
    path(reference_genome)
    path(reference_genome_index)
    path(baitset_bed)
    path(refflat_file)
    path(access_bed)
    path(targets_bed)
    path(antitargets_bed)
    path(baitset_genes_file)
    path(exclude_samples)

    output:
    tuple val(meta), path("parameters.json"), emit: parameters

    script:
    def exclude_arg = exclude_samples.name != 'NO_FILE' ? "-e ${exclude_samples}" : ''
    def contig_args = params.unplaced_contig_prefixes.collect { prefix -> "\"${prefix}\"" }.join(' ')
    """
    fur_cnvkit cnvkit_static_files \\
        -b ${bams} \\
        ${exclude_arg} \\
        -f ${reference_genome} \\
        -t ${baitset_bed} \\
        -u ${contig_args} \\
        -r ${refflat_file} \\
        -m ${metadata_manifest} \\
        -o . \\
        --metadata-sample-id-column "${params.col_sample_id}" \\
        --metadata-tumour-normal-column "${params.col_TN}" \\
        --metadata-sex-column "${params.col_sex}" \\
        --parameter-file-only \\
        --existing-access-bed ${access_bed} \\
        --existing-targets-bed ${targets_bed} \\
        --existing-antitargets-bed ${antitargets_bed} \\
        --existing-baitset-genes-file ${baitset_genes_file} \\
        --verbose DEBUG

    # Absolute paths in the parameter file point into this task directory, which no
    # longer exists once CNVKIT_CALL_COPY_NUMBER reads the file.
    rebase_parameter_paths.py --input parameters.json --output parameters.rebased.json
    mv parameters.rebased.json parameters.json
    """

    stub:
    def bam_names = bams.collect { bam -> "\"${bam.name}\"" }.join(', ')
    """
    cat <<JSON > parameters.json
    {
        "all_bams": [${bam_names}],
        "tumour_bams": [${bam_names}],
        "normal_bams": [],
        "reference_fasta": "${reference_genome.name}",
        "unplaced_contig_prefixes": [],
        "baitset_bed": "${baitset_bed.name}",
        "refflat_file": "${refflat_file.name}",
        "sample_metadata_xlsx": "${metadata_manifest.name}",
        "access_bed": "${access_bed.name}",
        "targets_bed": "${targets_bed.name}",
        "antitargets_bed": "${antitargets_bed.name}",
        "baitset_genes_file": "${baitset_genes_file.name}"
    }
    JSON
    """
}

process CNVKIT_POOLED_REFERENCE {
    label 'process_extralong'
    label 'process_high_memory'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    publishDir path: "${params.outdir}/pooled_reference",
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(parameters), path(bams), path(bais)
    path(metadata_manifest)
    path(reference_genome)
    path(reference_genome_index)
    path(baitset_bed)
    path(static_files)

    output:
    tuple val(meta), path("male/male.reference.cnn"), emit: male_reference
    tuple val(meta), path("female/female.reference.cnn"), emit: female_reference
    tuple val(meta), path("*/samples_used_in_reference.txt"), emit: samples_used
    tuple val(meta), path("*/normal_vs_normal"), emit: normal_vs_normal, optional: true

    script:
    """
    fur_cnvkit copy_number_reference \\
        -p ${parameters} \\
        -o . \\
        --max_cpus ${task.cpus}
    """

    stub:
    """
    mkdir -p male female
    echo stub > male/male.reference.cnn
    echo stub > female/female.reference.cnn
    echo stub > male/samples_used_in_reference.txt
    echo stub > female/samples_used_in_reference.txt
    """
}

process CNVKIT_CALL_COPY_NUMBER {
    label 'process_high'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    publishDir path: "${params.outdir}/cnvkit_cn_calling",
               mode: "${params.publish_dir_mode}",
               overwrite: true

    input:
    tuple val(meta), path(parameters), path(bams), path(bais)
    path(metadata_manifest)
    path(baitset_genes_file)
    path(male_reference)
    path(female_reference)

    output:
    tuple val(meta), path("*/*.call.median_centred.weight_ge_*.cns"), emit: cns
    tuple val(meta), path("*/*.metrics.tsv"), emit: metrics
    tuple val(meta), path("*/*.cnr"), emit: cnr
    tuple val(meta), path("*/*.{png,pdf}"), emit: qc_plots, optional: true

    script:
    // The thresholds are passed with `=` and no space: they start with a minus sign, so
    // `-t -0.737,...` would be parsed as another option by argparse.
    """
    fur_cnvkit cnvkit_cn_calling_pipeline \\
        -p ${parameters} \\
        -m ${male_reference} \\
        -f ${female_reference} \\
        -o . \\
        --weight_filter_threshold ${params.weight_filter_threshold} \\
        --call-thresholds="${params.call_thresholds}" \\
        --skip-genemetrics \\
        --verbose DEBUG
    """

    stub:
    // "default" is the study id fur-cnvkit derives from a TSV manifest; a real run reads
    // it from the manifest, and the output globs above match whatever it turns out to be.
    def thresholds = params.call_thresholds.replace(',', '_')
    def weight = params.weight_filter_threshold
    """
    mkdir -p default
    for bam in ${bams}; do
        sample=\$(echo \${bam} | cut -d. -f1)
        prefix="default/\${sample}.sample.dupmarked"
        echo stub > "\${prefix}.call.thresholds_${thresholds}.call.median_centred.weight_ge_${weight}.0.cns"
        echo stub > "\${prefix}.cnr"
    done
    echo stub > default/default.metrics.tsv
    """
}

process CNVKIT_RECALL {
    label 'process_low'
    container "gitlab-registry.internal.sanger.ac.uk/dermatlas/fur/fur_cnvkit:0.7.0"
    // The re-called file keeps the name of the file it was called from whenever the
    // thresholds are unchanged, so it is written into a sub-directory to avoid colliding
    // with its own input. That sub-directory is an artefact of the task, not of the
    // results, and is stripped on the way out.
    publishDir path: "${params.outdir}/cnvkit_cn_recall",
               mode: "${params.publish_dir_mode}",
               overwrite: true,
               saveAs: { filename -> filename.replaceFirst('^re-call/', '') }

    input:
    tuple val(meta), path(cns)
    val(recall_thresholds)

    output:
    tuple val(meta), path("re-call/*.cns"), emit: cns

    script:
    // cnvkit is re-run purely to rewrite the `cn` column: `--method none` leaves the
    // log2 ratios untouched, so the segments themselves are unchanged. `-y` tells cnvkit
    // the sample has one copy of each sex chromosome; without it a male sample's normal
    // single-copy X and Y are called as losses.
    def male_flag = meta.sex == 'M' ? '-y' : ''
    def default_underscored = params.call_thresholds.replace(',', '_')
    def recall_underscored = recall_thresholds.replace(',', '_')
    def output_name = cns.name.replace("thresholds_${default_underscored}", "thresholds_${recall_underscored}")
    """
    mkdir -p re-call
    cnvkit.py call ${cns} \\
        --method none \\
        -m threshold \\
        -t=${recall_thresholds} \\
        ${male_flag} \\
        -o re-call/${output_name}
    """

    stub:
    def default_underscored = params.call_thresholds.replace(',', '_')
    def recall_underscored = recall_thresholds.replace(',', '_')
    def output_name = cns.name.replace("thresholds_${default_underscored}", "thresholds_${recall_underscored}")
    """
    mkdir -p re-call
    echo stub > re-call/${output_name}
    """
}
