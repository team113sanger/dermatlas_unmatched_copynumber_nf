include { FILTER_SEGMENTS } from '../../modules/filter_segments.nf'
include { PENETRANCE_PLOT } from '../../modules/penetrance.nf'
include { RUN_GISTIC2 } from '../../modules/gistic2.nf'

// Everything that happens per subcohort: segment filtering, the GISTIC2 input .seg file,
// GISTIC2 itself, and the penetrance plot for the segments that went into it.
//
// Each subcohort is analysed twice, as a 'filtered' and an 'unfiltered' arm. The
// unfiltered arm is not a fallback - it is the control. Filtering drops segments that
// overlap difficult regions or are too short, and the only way to see whether it also
// dropped something real (a focal MYC amplification, a TP53 deletion) is to score both
// and compare.
workflow ANALYSE_SUBCOHORT {
    take:
    subcohorts           // tuple(meta, include_samples, exclude_samples) per subcohort
    recalled_cns         // value channel: list of re-called .cns files for the cohort
    difficult_regions
    gistic_refgene_file
    thresholds           // parsed threshold map (see modules/utils.nf)

    main:

    // One arm per (subcohort, analysis_type). The thresholds ride in the meta map because
    // both GISTIC2 (as amplitude cut-offs) and the penetrance plots (as log2 ratios) are
    // keyed on the same two numbers.
    arms = subcohorts
        .combine(Channel.of('filtered', 'unfiltered'))
        .map { meta, include_samples, exclude_samples, analysis_type ->
            tuple(
                meta + [
                    "analysis_type": analysis_type,
                    "label": analysis_type,
                    "plot_dir": "${meta.cohort_id}/${analysis_type}",
                    "gain": thresholds.gain,
                    "loss": thresholds.loss,
                    "loss_amplitude": thresholds.loss_amplitude
                ],
                include_samples,
                exclude_samples
            )
        }

    // The list of .cns files is wrapped so that combine() treats it as one element rather
    // than spreading each file into its own tuple slot.
    filter_input = arms
        .combine(recalled_cns.map { cns_files -> [cns_files] })
        .map { meta, include_samples, exclude_samples, cns_files ->
            tuple(meta, cns_files, include_samples, exclude_samples)
        }

    FILTER_SEGMENTS(filter_input, difficult_regions)

    // A subcohort can reach this point with samples listed and still have none left to
    // score: a small one (related_tumours, typically) can lose every sample to exclusion or
    // hypersegmentation. GISTIC2 cannot score an empty .seg file, so such an arm stops here
    // with a warning - its .seg and sample list are still published - rather than failing
    // the run for the subcohorts that have samples.
    scored_samples = FILTER_SEGMENTS.out.samples
        .filter { meta, samples ->
            def has_samples = samples.readLines().any { line -> line.trim() }
            if (!has_samples) {
                log.warn("Skipping GISTIC2 and penetrance plot for subcohort '${meta.cohort_id}' " +
                         "(${meta.analysis_type}): no samples left after exclusion and filtering.")
            }
            return has_samples
        }

    RUN_GISTIC2(
        FILTER_SEGMENTS.out.seg
            .join(scored_samples)
            .map { meta, seg, _samples -> tuple(meta, seg) },
        gistic_refgene_file
    )

    // The sample list FILTER_SEGMENTS wrote is the definitive membership of the .seg file
    // it exported - the subcohort's samples, minus the excluded and hypersegmented ones -
    // so the plot is drawn over exactly the samples GISTIC2 scored.
    //
    // Only the filtered arm rewrites its .cns files; the unfiltered arm exported the
    // re-called ones untouched, so `remainder: true` keeps it in the stream with a null
    // where its .cns files would be, and it is plotted from the re-called files instead.
    // An arm skipped above has .cns files but no scored samples; the remainder join would
    // emit it with a null sample list, so it is dropped here.
    segments_to_plot = scored_samples
        .join(FILTER_SEGMENTS.out.cns, remainder: true)
        .filter { _meta, samples, _cns -> samples != null }
        .branch { meta, samples, cns ->
            filtered: cns != null
                return tuple(meta, cns, samples)
            unfiltered: true
                return tuple(meta, samples)
        }

    unfiltered_to_plot = segments_to_plot.unfiltered
        .combine(recalled_cns.map { cns_files -> [cns_files] })
        .map { meta, samples, cns_files -> tuple(meta, cns_files, samples) }

    PENETRANCE_PLOT(
        segments_to_plot.filtered
            .mix(unfiltered_to_plot)
            .map { meta, cns, samples -> tuple(meta, cns, samples, file("${projectDir}/assets/NO_FILE")) }
    )

    emit:
    seg            = FILTER_SEGMENTS.out.seg
    filtered_cns   = FILTER_SEGMENTS.out.cns
    removed_cns    = FILTER_SEGMENTS.out.removed
    gistic_lesions = RUN_GISTIC2.out.lesions
    gistic_scores  = RUN_GISTIC2.out.scores
    plots          = PENETRANCE_PLOT.out.plot
}
