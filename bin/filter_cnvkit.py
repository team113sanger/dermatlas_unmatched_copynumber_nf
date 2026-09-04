#!/usr/bin/env python3
import pandas as pd
import argparse
import logging
import os
import inspect
from pathlib import Path
from skgenome import tabio
from cnvlib.export import export_seg

# ----- Configure logging -----

logging.basicConfig(
    level=logging.INFO
)

# ----- Define functions -----

def get_bed_intervals(bed_file):
    # Read in BED regions for filtering
    bed_intervals = {}
    with open(bed_file) as bed:
        for line in bed:
            if (
                line.strip() == ""
                or line.startswith("#")
                or line.startswith("chrom")
            ):
                continue
            chrom, start, end, *rest = line.rstrip().split("\t")
            start, end = int(start), int(end)

            bed_intervals.setdefault(chrom, []).append((start, end))

    return bed_intervals


def filter_cns(cns_files, bed_intervals=None, threshold=0.5, weight=0, remove_seg=0, ext="filtered", out_dir=None):
    new_cns_files = []
    for file in cns_files:
        output_lines = []
        removed_lines = []

        p = Path(file)
        ext = ext.lstrip(".")
        removed_ext = "removed"

        new_name = f"{p.stem}.{ext}{p.suffix}"
        removed_name = f"{p.stem}.{removed_ext}{p.suffix}"

        out_path = Path(out_dir) if out_dir else p.parent
        out_path.mkdir(parents=True, exist_ok=True)

        output_file = out_path / new_name
        removed_file = out_path / removed_name

        with open(file) as cns:
            logging.info(f"PROCESSING FILE {file}")
            logging.info(f"Filtering by weight < {weight}")
            if bed_intervals is not None:
                logging.info(f"Overlapping with intervals, threshold is {threshold}")
            else:
                logging.info(f"Not overlapping with BED intervals")

            header = cns.readline().strip().split("\t")
            col_index = {name: idx for idx, name in enumerate(header)}

            chrom_idx = col_index['chromosome']
            start_idx = col_index['start']
            end_idx = col_index['end']
            weight_idx = col_index['weight']

            output_lines.append("\t".join(header) + "\n")
            removed_lines.append("\t".join(header) + "\n")

            for line in cns:
                remove = False
                parts = line.rstrip().split("\t").copy()

                seg_start, seg_end = int(parts[start_idx]), int(parts[end_idx])
                seg_len = seg_end - seg_start

                if bed_intervals is not None:
                    chrom = parts[chrom_idx]

                    if chrom in bed_intervals:
                        tot_overlap = 0
                        for bstart, bend in bed_intervals[chrom]:
                            if bend < seg_start:
                                continue
                            if bstart > seg_end:
                                break

                            # Compute overlap length
                            overlap = max(0, min(seg_end, bend) - max(seg_start, bstart))
                            tot_overlap += overlap

                        # Remove if overlap ≥ X% of the segment
                        if tot_overlap >= threshold * seg_len:
                            remove = True
                            print(f"Removing by overlap: {tot_overlap} {seg_len} {line}")
#                            break
                        else:
                            print(f"Keeping: {tot_overlap} {seg_len} {line}")

                if weight > 0:
                    seg_weight = float(parts[weight_idx])
                    if seg_weight < weight:
                        remove = True
                        print(f"Removing by weight filter {line}")
                if remove_seg > 0:
                    if seg_len < remove_seg:
                        remove = True
                        print(f"Removing by segment sie filter {line}")

                if remove:
                    removed_lines.append(line)
                else:
                    output_lines.append(line)

        # Write output files
        with open(output_file, "w") as out:
            logging.info(f"Writing filtered lines to: {output_file}")
            #for l in output_lines:
            #    out.write(l)
            out.writelines(output_lines)

        with open(removed_file, "w") as out:
            logging.info(f"Writing removed lines to: {removed_file}")
            out.writelines(removed_lines)

        new_cns_files.append(output_file)

    return new_cns_files


def export_gistic(cns_files, output, samples=None, strip_filenames=False):
    df = export_seg(cns_files)

    # Filter table to keep only samples listed
    if samples is not None:
        logging.info(f"Samples used for GISTIC *.seg file: {samples}")
        # Filter dataframe
        df['sample_id'] = df['ID'].str.split('.').str[0]
        filtered_df = df[df['sample_id'].isin(samples)].copy()
        
        # Strip filenames if requested
        if strip_filenames:
            logging.info(f"Stripping file names in GISTIC *.seg file")
            filtered_df['ID'] = filtered_df['sample_id']
        filtered_df = filtered_df.drop('sample_id', axis=1)
    else:
        logging.info(f"Using all samples from *cns files for GISTIC *.seg")

    filtered_df.to_csv(output, sep='\t', index=False)
    logging.info(f"Filtered segments saved to: {output}")

    #return seg_table


def main():

    parser = argparse.ArgumentParser(
        usage="python %(prog)s --input *median-centred.cns --samples samples.list [options]",
        description=(
            "Filter CNVkit *cns files and create GISTIC2 .seg files.\n"
            "\n"
            "Filter *cns files by one or all of: segment weight, quality metrics\n"
            "(mean segement size), overlap with BED file regions."
        ),
        epilog="""Examples:
        (1) Filter the *cns files by segment weight, quality metrics and BED overlap.
            Overlap must be at least 25%% of segment size.

        python %(prog)s \\
            --input cnvkit_cn_calling/Sheet1/*.call.median_centred.cns \\
            --samples samples_list.txt \\
            --filter-cns \\
            --weight 30 \\
            --suffix filtered \\
            --bed-file mask_regions.bed \\
            --bed-threshold 0.25 \\
            --quality cnvkit_cn_calling/Sheet1/Sheet1.metrics.tsv \\
            --sample-list-out samples_to_use.txt

        (2) Filter the *cns files by segment weight and mean seg. size (quality metrics file)
            and create a GISTIC2 *seg file from the filtered *cns files

        python %(prog)s \\
            --input cnvkit_cn_calling/Sheet1/*.call.median_centred.cns \\
            --filter-cns \\
            --weight 30 \\
            --quality cnvkit_cn_calling/Sheet1/Sheet1.metrics.tsv \\
            --suffix filtered \\
            --samples samples_list.txt \\
            --sample-list-out samples_to_use.txt \\
            --gistic-seg gistic_inputs/gistic_filtered.seg

        (3) Do not filter *cns files, but filter samples using mean seg. size
            (quality metrics file) and create a GISTIC2 *.seg file.

        python %(prog)s \\
            --input cnvkit_cn_calling/Sheet1/*.call.median_centred.cns \\
            --samples samples_list.txt \\
            --quality cnvkit_cn_calling/Sheet1/Sheet1.metrics.tsv \\
            --suffix filtered \\
            --sample-list-out samples_to_use.txt \\
            --gistic-seg gistic_inputs/gistic_filtered.seg

    """,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Main input
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i', nargs="+", required=True, help='Input segments files (.cns)')
    required.add_argument('--samples', '-s', required=True, help='Samples to include file')
    
    # cns file filtering options
    filter_group = parser.add_argument_group('Filtering options for *.cns files')
    
    filter_group.add_argument('--filter-cns', '-f', action="store_true", help='Filter cns file segements. Use with --bed and/or --weight. Optional')

    # Filter by overlap with regions in BED file
    filter_group.add_argument('--bed-file', '-b', required=False, help='Optional BED file of regions to exclude')
    filter_group.add_argument(
        '--bed-threshold', '-t',
        type=float,
        required=False,
        default = 0.5,
        help='Minimumn fraction segment overlap with BED file region for exclusion from cns file (default = 0.5)'
    )

    # Filter by segment weight value
    filter_group.add_argument(
        '--weight', '-w',
        type=float,
        required=False,
        default = 0,
        help='A cns weight threshold for cns filtering. Use with --filter-cns (default: 0)'
    )

    # Remove segments by segment size
    filter_group.add_argument(
        '--remove-seg',
        type=int,
        required=False,
        default = 0,
        help='A minimumn threshold for segment size. Use with --filter-cns (default: 0)'
    )

    # Filter samples by pre-calculated average segment size in quality metrics file
    filter_group.add_argument('--quality', '-q', required=False, help='Quality metrics file')
    filter_group.add_argument('--min-seg-size', type=int, default=4000000, 
                        help='Minimum average segment size per sample (default: 4000000)')

    # Output file options
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument('--exclude-samples', '-e', required=False, help='Samples to exclude from GISTIC2 seg files')
    output_group.add_argument(
        '--suffix',
        required=False,
        default='filtered',
        help='Suffix to add to cns/seg file name before file extention'
    )
    output_group.add_argument('--sample-list-out', '-sl', required=False, help='Output samples in the filtered set')
    output_group.add_argument('--output-dir', required=False, help='Path to output directory for filtered cns files')

    # GISTIC export options
    output_group.add_argument('--gistic-seg', required=False, help='Create a gistic2 .seg file from *cns inputs')
    output_group.add_argument('--strip-filenames', action='store_true',
                        help='Replace ID with sample_id in GISTIC2 .seg file (strip file extensions)'),

    args = parser.parse_args()

    # Check input parameters

    print("\nSELECTED OPTIONS:")
    maxlen = max(len(k) for k in vars(args))

    for key, value in vars(args).items():
        print(f"  {key.ljust(maxlen)} : {value}")

    print("\n")

    if (args.quality or args.gistic_seg) and not args.sample_list_out:
        parser.error("--sample-list-out must be provided if --quality or --gistic-seg is provided")

    if args.filter_cns and not (args.bed_file or args.weight or args.remove_seg):
        parser.error("If --filter-cns is provided, use one or all of --bed-file or --weight or --remove-seg")

    if (args.bed_file or args.weight or args.remove_seg) and not args.filter_cns:
        parser.error("Use --filter-cns if filtering with --bed-file and/or --weight")

    # Read the segments (cns) files

    used_cns_files = args.input

    if args.filter_cns:
        # Get BED intervals to remove from cns
        if args.bed_file is not None:
            logging.info(f"Reading BED file intervals from {args.bed_file}")
            intervals = get_bed_intervals(args.bed_file) 
        else:
            logging.info(f"No BED file provided - not filtering by BED regions")
            intervals = None
            args.bed_threshold = "NA"
        if args.weight == 0:
            logging.info(f"Not filtering by segment weight")
        else:
            logging.info("Filtering: keeping segments with weight >= {args.weight}")

        # Filter cns files
        logging.info(f"Beginning *cns filtering: bed_threshold={args.bed_threshold} weight={args.weight} file suffix={args.suffix}")
        used_cns_files = filter_cns(args.input, intervals, args.bed_threshold, args.weight, args.remove_seg, args.suffix, args.output_dir)

    # Read samples to check/include
    logging.info(f"Reading sample file {args.samples}")
#    samples_to_include = Path(args.samples).read_text().splitlines()
    samples_to_include = [
        line.split()[0]
        for line in Path(args.samples).read_text().splitlines()
        if line.strip()
    ]

    samples_to_exclude = []

    if args.exclude_samples is not None:
        logging.info(f"Reading sample exclude file {args.exclude_samples}")
        samples_to_exclude = [
            line.split()[0]
            for line in Path(args.exclude_samples).read_text().splitlines()
            if line.strip()
        ]

    # Read quality metrics and create a list of samples that pass
    logging.info(f"Number of samples to check (sample list): {len(samples_to_include)}")
    logging.info(f"Number of cns files: {len(used_cns_files)}")
    print(used_cns_files)
    cns_samples = [Path(fname).name.split('.', 1)[0] for fname in used_cns_files]
    #print(cns_samples)
    cns_samples_to_include = [s for s in cns_samples if s in samples_to_include]    
    #print(cns_samples_to_include)
    logging.info(f"Number of sample cns files in sample list: {len(cns_samples_to_include)}")

    if args.quality is not None:
        logging.info(f"Applying sample filtering using {args.quality} with min_seg_size={args.min_seg_size}")
        quality_metrics = pd.read_csv(args.quality, sep='\t')
        high_quality_samples = quality_metrics[quality_metrics['AvgSegSize'] > args.min_seg_size]
        low_quality_samples = quality_metrics[quality_metrics['AvgSegSize'] <= args.min_seg_size]
        high_quality_sample_ids = high_quality_samples['Sample'].tolist()
        samples_to_exclude.extend(low_quality_samples['Sample'].tolist())

        # Filter samples by quality
        samples_filtered = [s for s in cns_samples_to_include if s in high_quality_sample_ids]
        logging.info(f"Hypersegmented samples removed: {len(cns_samples_to_include) - len(samples_filtered)}")

        cns_samples_to_include = samples_filtered

        # File for hypersegmented samples
        logging.info(f"Total unique samples removed (excluded+hypersegmented): {len(samples_to_exclude)}")
    else:
        logging.info(f"Not checking quality metrics for hypersegmented samples")
        logging.info(f"Excluded samples provided by user (with or without cns files): {len(samples_to_exclude)}")
 
    samples_to_exclude = list(set(samples_to_exclude))        
    cns_samples_to_include = [s for s in cns_samples_to_include if s not in samples_to_exclude]    
    outpath = Path(args.sample_list_out)
    exclude_file = outpath.with_name(f"{outpath.stem}.excluded{outpath.suffix}")
    pd.Series(samples_to_exclude).to_csv(
        exclude_file, sep='\t', index=False, header=False)
    logging.info(f"Excluded sample list saved to: {exclude_file}")
    logging.info(f"Final included sample count: {len(cns_samples_to_include)}")

    # Save either full list or filtered list

    pd.Series(cns_samples_to_include).to_csv(
        args.sample_list_out, sep='\t', index=False, header=False)
    logging.info(f"Final sample list saved to: {args.sample_list_out}")

    # Create gistic2 inputs
    if args.gistic_seg is not None:
        if args.filter_cns:
            logging.info(f"Using *cns files with additional filtering for GISTIC *.seg")
        else:
            logging.info(f"Using *cns files without additional for GISTIC *.seg")

        # Create GISTIC2 inputs
        export_gistic(used_cns_files, args.gistic_seg, samples=cns_samples_to_include, strip_filenames=args.strip_filenames)


# ----- Entry point ----- 
    
if __name__ == "__main__":
    main()
