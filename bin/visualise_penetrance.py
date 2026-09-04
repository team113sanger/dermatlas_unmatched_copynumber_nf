#!/usr/bin/env python3
"""
Copy Number Penetrance Analysis CLI

A command-line tool for generating genome-wide copy number penetrance plots
from CNVkit segmentation data.
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path
from typing import List, Tuple, Optional


def generate_penetrance_df(
    cns_files: List[Path], 
    bin_size: int = 1_000_000, 
    gain_threshold: float = 0.32, 
    loss_threshold: float = -0.4
) -> pd.DataFrame:
    """
    Generate penetrance dataframe from CNS files.
    
    Args:
        cns_files: List of paths to .cns files
        bin_size: Size of genomic bins in base pairs
        gain_threshold: Log2 ratio threshold for gains
        loss_threshold: Log2 ratio threshold for losses
        
    Returns:
        DataFrame with penetrance data
    """
    all_bins = []

    num_samples = len(cns_files)
 
    for file in cns_files:
        try:
            sample_data = pd.read_csv(file, sep='\t')
            sample_id = file.name.split('.')[0]
            
            for index, row in sample_data.iterrows():
                chrom = row['chromosome']
                start = row['start']
                end = row['end']
                log2_ratio = row['log2']
                
                bin_start = (start // bin_size) * bin_size
                bin_end = ((end // bin_size) + 1) * bin_size
                
                for bin_pos in range(bin_start, bin_end, bin_size):
                    overlap_start = max(start, bin_pos)
                    overlap_end = min(end, bin_pos + bin_size)
                    
                    if overlap_start < overlap_end:
                        overlap_fraction = (overlap_end - overlap_start) / (end - start)
                        
                        all_bins.append({
                            'Chromosome': chrom,
                            'Bin_Start': bin_pos,
                            'Bin_End': bin_pos + bin_size,
                            'Sample_ID': sample_id,
                            'Gain': 1 if log2_ratio > gain_threshold else 0,
                            'Loss': 1 if log2_ratio < loss_threshold else 0,
                            'Overlap_Fraction': overlap_fraction
                        })
        except Exception as e:
            print(f"Warning: Error processing {file}: {e}", file=sys.stderr)
            continue
    
    if not all_bins:
        raise ValueError("No valid data found in input files")
    
    bins_df = pd.DataFrame(all_bins)
    
    # Calculate penetrance for bins with data
    penetrance_df = bins_df.groupby(['Chromosome', 'Bin_Start', 'Bin_End']).agg(
        Gain_Penetrance=('Gain', 'sum'),
        Loss_Penetrance=('Loss', 'sum'),
        #DEBUG
        Gain_Samples=('Sample_ID',
                      lambda x: samples_to_string(
                          x[bins_df.loc[x.index, 'Gain'] > 0]
                      )),
        Loss_Samples=('Sample_ID',
                      lambda x: samples_to_string(
                          x[bins_df.loc[x.index, 'Loss'] > 0]
                      ))
    ).reset_index()

    penetrance_df['Gain_Penetrance'] = penetrance_df['Gain_Penetrance'] * 100 / num_samples
    penetrance_df['Loss_Penetrance'] = penetrance_df['Loss_Penetrance'] * 100 / num_samples

    # Create complete grid of bins for each chromosome (from 0 to max position)
    complete_bins = []
    for chrom in penetrance_df['Chromosome'].unique():
        chrom_data = penetrance_df[penetrance_df['Chromosome'] == chrom]
        max_pos = chrom_data['Bin_End'].max()
        
        # Create all bins from 0 to max_pos
        for bin_start in range(0, max_pos, bin_size):
            complete_bins.append({
                'Chromosome': chrom,
                'Bin_Start': bin_start,
                'Bin_End': bin_start + bin_size
            })
    
    complete_df = pd.DataFrame(complete_bins)
    
    # Merge with actual data, filling missing bins with 0
    penetrance_df = complete_df.merge(
        penetrance_df,
        on=['Chromosome', 'Bin_Start', 'Bin_End'],
        how='left'
    )
    
    # Fill NaN values with 0 (bins with no data)
    penetrance_df['Gain_Penetrance'] = penetrance_df['Gain_Penetrance'].fillna(0)
    penetrance_df['Loss_Penetrance'] = penetrance_df['Loss_Penetrance'].fillna(0)
    
    #penetrance_df['Gain_Penetrance'] *= 100
    #penetrance_df['Loss_Penetrance'] *= 100
    
    return penetrance_df

def plot_penetrance(
    data: pd.DataFrame,
    n_samples: int,
    title: str = 'Genome-Wide DNA Copy Number Penetrance',
    bin_size: int = 1_000_000,
    gain_threshold: float = 0.32,
    loss_threshold: float = -0.4,
    y_range: Tuple[int, int] = (-50, 50),
    y_tick_step: int = 10,
    output_file: Optional[str] = None,
    exclude_chromosomes: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (30, 6),
    dpi: int = 300,
    fontsize: int = 18,
    label_fontsize: int = 12
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot genome-wide copy number penetrance.
    
    Args:
        data: Penetrance dataframe
        title: Plot title
        bin_size: Genomic bin size
        gain_threshold: Gain threshold for legend
        loss_threshold: Loss threshold for legend
        y_range: Y-axis range as (min, max)
        y_tick_step: Y-axis tick step size
        output_file: Output file path (if None, shows plot)
        exclude_chromosomes: List of chromosomes to exclude
        figsize: Figure size as (width, height)
        dpi: Output DPI for saved plots
        
    Returns:
        Tuple of (figure, axes) objects
    """
    if exclude_chromosomes is None:
        exclude_chromosomes = []
    
    exclude_chromosomes = [str(chrom).replace('chr', '') for chrom in exclude_chromosomes]
    data['Chromosome'] = data['Chromosome'].astype(str)
    data = data[~data['Chromosome'].str.replace('chr', '').isin(exclude_chromosomes)]
    
    unique_chromosomes = list(data['Chromosome'].unique())
    
    def chr_sort_key(chrom):
        chrom_clean = chrom.replace('chr', '')
        if chrom_clean.isdigit():
            return (0, int(chrom_clean))
        elif chrom_clean == 'X':
            return (1, 23)
        elif chrom_clean == 'Y':
            return (1, 24)
        elif chrom_clean in ['MT', 'M']:
            return (1, 25)
        else:
            return (2, chrom_clean)
    
    chromosomes = sorted(unique_chromosomes, key=chr_sort_key)
    print(f"Processing chromosomes: {chromosomes}")
    
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    
    x_pos = 0
    x_ticks = []
    x_labels = []
    chromosome_boundaries = []
    
    gain_label = f'Gain (log2 > {gain_threshold})'
    loss_label = f'Loss (log2 < {loss_threshold})'
    
    for chrom in chromosomes:
        chrom_data = data[data['Chromosome'] == chrom]
        
        if chrom_data.empty:
            continue
        
        x = chrom_data['Bin_Start'] + x_pos - chrom_data['Bin_Start'].min()
        
        ax.fill_between(
            x, 0, chrom_data['Gain_Penetrance'],
            color='red', alpha=1,
            label=gain_label if chrom == chromosomes[0] else ""
        )
        ax.fill_between(
            x, 0, -chrom_data['Loss_Penetrance'],
            color='blue', alpha=1,
            label=loss_label if chrom == chromosomes[0] else ""
        )
        
        length = chrom_data['Bin_End'].max() - chrom_data['Bin_Start'].min() + bin_size
        x_pos += length
        
        x_ticks.append(x_pos - length / 2)
        x_labels.append(chrom)
        chromosome_boundaries.append(x_pos)
    
    ax.set_xlim(0, chromosome_boundaries[-1])
    
    y_min, y_max = y_range
    ax.set_ylim(y_min, y_max)
    ax.set_yticks(range(y_min, y_max + 1, y_tick_step))
    ax.tick_params(axis='y', labelsize=label_fontsize)
    
    ax.set_title(f"{title} (n = {n_samples})", fontsize=fontsize + 2)
    ax.set_xlabel('Genomic Position (by Chromosome)', fontsize=fontsize)
    ax.set_ylabel('Penetrance (%)', fontsize=fontsize)
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=label_fontsize)
    for b in chromosome_boundaries + [0]:
        ax.axvline(x=b, color='black', linestyle='-', linewidth=0.8)
    ax.grid(True, which='both', axis='y', linestyle='--', linewidth=0.5)
    
#    ax.legend(
#        loc='center left',
#        bbox_to_anchor=(1.02, 0.5),
#        borderaxespad=0,
#        fontsize=fontsize
#    )
    
    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.2),
        ncol=2,
        borderaxespad=0,
        fontsize=fontsize
    )
    
    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    else:
        plt.show()
    
    return fig, ax


def load_sample_list(file_path: Path) -> List[str]:
    """Load sample list from file."""
    if not file_path.exists():
        raise FileNotFoundError(f"Sample list file not found: {file_path}")
    return file_path.read_text().strip().splitlines()


def find_cns_files(base_dir: Path, pattern: str = "*.median_centred.cns") -> List[Path]:
    """Find CNS files in directory."""
    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory not found: {base_dir}")
    
    files = list(base_dir.glob(pattern))
    if not files:
        # Try alternative patterns
        alt_patterns = ["*.median_centred.cns"]
        for alt_pattern in alt_patterns:
            files = list(base_dir.glob(alt_pattern))
            if files:
                print(f"Found files with pattern: {alt_pattern}")
                break
    
    return files


def samples_to_string(series):
    return ';'.join(sorted(pd.unique(series)))


def main():
    parser = argparse.ArgumentParser(
        description="Generate genome-wide copy number penetrance plots from CNVkit data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python visualise_penetrance.py -d /path/to/cns/files -o penetrance_plot.pdf
  
  # With sample filtering
  python visualise_penetrance.py -d /path/to/cns/files -s samples.txt -e exclude.txt -o plot.pdf
  
  # Custom thresholds and parameters
  python visualise_penetrance.py -d /path/to/cns/files -o plot.pdf --gain-threshold 0.5 --loss-threshold -0.5 --bin-size 500000
        """
    )
    
    # Required arguments
    parser.add_argument('-d', '--data-dir', type=Path, required=True,
                       help='Directory containing CNS files')
    parser.add_argument('-o', '--output', type=Path, required=True,
                       help='Output plot file (PDF, PNG, SVG, etc.)')
    
    # Optional filtering
    parser.add_argument('-s', '--samples', type=Path,
                       help='File containing samples to include (one per line)')
    parser.add_argument('-e', '--exclude-samples', type=Path,
                       help='File containing samples to exclude (one per line)')
    parser.add_argument('--pattern', default='*median_centred.cns',
                       help='File pattern for CNS files (default: *.cns)')
    
    # Analysis parameters
    parser.add_argument('--bin-size', type=int, default=1_000_000,
                       help='Genomic bin size in bp (default: 1000000)')
    parser.add_argument('--gain-threshold', type=float, default=0.32,
                       help='Log2 ratio threshold for gains (default: 0.32)')
    parser.add_argument('--loss-threshold', type=float, default=-0.4,
                       help='Log2 ratio threshold for losses (default: -0.4)')
    
    # Plot parameters
    parser.add_argument('--title', default='Genome-Wide Copy Number Penetrance',
                       help='Plot title')
    parser.add_argument('--exclude-chromosomes', nargs='*', default=['X', 'Y'],
                       help='Chromosomes to exclude from plot (default: X Y)')
    parser.add_argument('--y-range', type=int, nargs=2, default=[-50, 50],
                       metavar=('MIN', 'MAX'),
                       help='Y-axis range (default: -50 50)')
    parser.add_argument('--y-tick-step', type=int, default=10,
                       help='Y-axis tick step size (default: 10)')
    parser.add_argument('--figsize', type=int, nargs=2, default=[30, 6],
                       metavar=('WIDTH', 'HEIGHT'),
                       help='Figure size in inches (default: 30 6)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Output DPI (default: 300)')
    parser.add_argument('--fontsize', type=float, default=18,
                       help='Axis/legend title font size; title is fontsize+2 (default: 18)')
    parser.add_argument('--label-fontsize', type=float, default=12,
                       help='Axis tick font size (default: 12)')
    
    # Output options
    parser.add_argument('--save-data', type=Path,
                       help='Save penetrance data to CSV file')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    try:
        # Find CNS files
        if args.verbose:
            print(f"Searching for files in: {args.data_dir}")
        cns_files = find_cns_files(args.data_dir, args.pattern)
        
        if not cns_files:
            print(f"Error: No CNS files found in {args.data_dir} with pattern {args.pattern}")
            sys.exit(1)
        
        if args.verbose:
            print(f"Found {len(cns_files)} CNS files")
            #print(f"{cns_files}")
            print("\n".join(map(str, cns_files)))
        
        # Load sample lists
        samples_to_include = None
        if args.samples:
            samples_to_include = load_sample_list(args.samples)
            if args.verbose:
                print(f"Loaded {len(samples_to_include)} samples to include")
        
        samples_to_exclude = []
        if args.exclude_samples:
            samples_to_exclude = load_sample_list(args.exclude_samples)
            if args.verbose:
                print(f"Loaded {len(samples_to_exclude)} samples to exclude")
        
        # Filter files
        filtered_files = []
        for file in cns_files:
            sample_id = file.name.split('.')[0]
            
            # Check inclusion list
            if samples_to_include and sample_id not in samples_to_include:
                continue
            
            # Check exclusion list
            if sample_id in samples_to_exclude:
                continue
            
            filtered_files.append(file)
        
        print(f"Processing {len(filtered_files)} valid samples...")
        
        if not filtered_files:
            print("Error: No valid samples found after filtering")
            sys.exit(1)
        
        # Generate penetrance data
        if args.verbose:
            print("Generating penetrance data...")
        
        penetrance_df = generate_penetrance_df(
            filtered_files,
            bin_size=args.bin_size,
            gain_threshold=args.gain_threshold,
            loss_threshold=args.loss_threshold
        )
        
        if args.verbose:
            print(f"Generated penetrance data with {len(penetrance_df)} bins")
        
        # Save data if requested
        if args.save_data:
            penetrance_df.to_csv(args.save_data, index=False)
            print(f"Penetrance data saved to: {args.save_data}")
        
        # Create output directory
        args.output.parent.mkdir(parents=True, exist_ok=True)
        
        # Generate plot
        if args.verbose:
            print("Generating plot...")
        
        fig, ax = plot_penetrance(
            penetrance_df,
            n_samples=len(filtered_files),
            title=args.title,
            bin_size=args.bin_size,
            gain_threshold=args.gain_threshold,
            loss_threshold=args.loss_threshold,
            y_range=tuple(args.y_range),
            y_tick_step=args.y_tick_step,
            output_file=str(args.output),
            exclude_chromosomes=args.exclude_chromosomes,
            figsize=tuple(args.figsize),
            dpi=args.dpi,
            fontsize=args.fontsize,
            label_fontsize=args.label_fontsize
        )

        print(f"Analysis complete! Plot saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
