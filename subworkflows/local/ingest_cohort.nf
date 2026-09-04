include { countDataRows; readHeader } from '../../modules/utils.nf'

// Turn the cohort's BAM directory and metadata manifest into the two channels the rest
// of the pipeline works from: one cohort-level bundle of BAMs, and a per-sample sex
// lookup.
//
// Every sample-level decision the original wrapper scripts made with `dir -1 bams | grep
// -v -F -f ...` happens here instead, before any compute is scheduled: the universe of
// analysable tumours restricts which BAMs are read at all, and the explicit exclude list
// removes named samples on top of that.
workflow INGEST_COHORT {
    take:
    bam_path            // glob matching indexed BAMs and their .bai
    metadata_manifest   // TSV manifest with sample id, sex and tumour/normal columns
    all_samples         // optional TSV listing every analysable tumour, or null
    excluded_samples    // list of sample ids to drop entirely (possibly empty)
    cohort_meta         // map carried alongside the cohort-level bundle

    main:

    // The glob must match the BAM and its index and nothing else: `size: 2` silently
    // drops any key that does not match exactly two files, so a stray .bam.bas sibling
    // turns a sample into no sample at all. Use "<dir>/**bam{,.bai}".
    bams = Channel.fromFilePairs(bam_path, size: 2, flat: true) { bam ->
            bam.name.tokenize('.').first()
        }
        .map { sample_id, bam, bai -> tuple(["sample_id": sample_id], bam, bai) }

    // Optional sample universe. Without it every BAM matched by the glob is analysed, so
    // a BAM left behind in the directory becomes a sample; with it, only the tumours the
    // study has agreed to analyse are read.
    universe = all_samples
        ? (file(all_samples, checkIfExists: true)
            .readLines()
            .collect { line -> line.trim().split('\t')[0].trim() }
            .findAll { id -> id } as Set)
        : null

    if (universe != null && universe.isEmpty()) {
        error("No sample ids read from --all_samples '${all_samples}'. Expected one sample id per line.")
    }

    excluded = excluded_samples as Set

    selected_bams = bams.filter { meta, bam, bai ->
        if (excluded.contains(meta.sample_id)) {
            return false
        }
        return universe == null || universe.contains(meta.sample_id)
    }

    // Warn about ids the study expects but the BAM directory does not hold: without this
    // a mistyped path or an incomplete transfer produces a smaller cohort that still
    // runs to completion.
    if (universe != null) {
        selected_bams
            .map { meta, bam, bai -> meta.sample_id }
            .collect()
            .subscribe { found ->
                def missing = (universe - excluded - (found as Set)).sort()
                if (missing) {
                    log.warn("${missing.size()} sample(s) listed in --all_samples have no BAM matching " +
                             "--bam_path and will not be analysed: ${missing.join(', ')}")
                }
            }
    }

    // One bundle per cohort: fur-cnvkit takes the whole cohort in a single invocation, so
    // the per-sample tuples are collected back into a pair of lists here.
    cohort_bams = selected_bams
        .map { meta, bam, bai -> tuple(bam, bai) }
        .collect(flat: false)
        .map { pairs ->
            tuple(cohort_meta, pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] })
        }

    // Sex per tumour sample, as read from the manifest. Only the tumours are kept: the
    // normals are in the manifest to build a pooled reference, and are never re-called.
    sample_sex = Channel.fromPath(metadata_manifest, checkIfExists: true)
        .splitCsv(sep: '\t', header: true)
        .filter { row -> row[params.col_TN] && row[params.col_TN] != 'N' }
        .map { row -> tuple(row[params.col_sample_id], row[params.col_sex]) }
        .filter { sample_id, sex ->
            if (!(sex in ['M', 'F'])) {
                // fur-cnvkit infers a sex for calling, but `cnvkit call` needs to be told
                // explicitly whether to expect one or two sex chromosomes, so a sample of
                // unknown sex cannot be re-called.
                log.warn("Sample '${sample_id}' has sex '${sex}' in the manifest; it will be called but not re-called.")
                return false
            }
            return true
        }

    emit:
    bams        = selected_bams
    cohort_bams = cohort_bams
    sample_sex  = sample_sex
}
