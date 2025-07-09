#!/usr/bin/env python3
"""
A Python console program to analyze BAM files.
This script processes BAM files,
-- generates plots of
the percentage of mapped calls by read position
  and
the distribution of PHRED scores by read position.
"""
import argparse
import os
import subprocess
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pysam

# noinspection PyUnresolvedReferences
import seaborn as sns


def wget_download_file_if_not_exists(url, path_local: "os.PathLike[str] | str"):
    """Download a file from a URL to a local path if it doesn’t exist.

    Args:
        url (str): The URL of the file to download.
        path_local (os.PathLike[str] | str): The local path to save the downloaded file.
    """
    if not os.path.exists(path_local):
        print(f'Downloading "{url}" to "{path_local}" ...')
        subprocess.run(["wget", str(url), "-O", str(path_local)], check=True)
    else:
        print(f'Using existing file "{path_local}" .')


def main(chromosome_number, base_url, bam_path, read_length):
    """
    Main function to process a BAM file and perform analyses.

    This function
    opens a BAM file,
    initializes data structures for
    mapped call counts and PHRED scores, and
    processes all reads aligned to chromosome 20.

    It calculates and plots
    the percentage of mapped calls by read position
      and
    the distribution of PHRED scores by read position.
    """

    url = f"{base_url}{bam_path.name}"

    wget_download_file_if_not_exists(url, bam_path)
    wget_download_file_if_not_exists(f"{url}.bai", f"{bam_path}.bai")

    bam = pysam.AlignmentFile(bam_path, 'rb')

    # Initialize data structures
    counts = [0] * read_length  # For mapped call counts per read position (0 - (read_length-1))
    phreds = defaultdict(list)  # For PHRED scores per read position
    total_reads = 0

    # Process all reads aligned to chosen chromosome in a single pass
    print(f"Processing {bam_path} on chromosome {chromosome_number}:")
    for rec in bam.fetch(f'{chromosome_number}'):
        total_reads += 1
        for i in range(rec.query_alignment_start, rec.query_alignment_end):
            counts[i] += 1
            phreds[i].append(rec.query_qualities[i])
        print(f"\r    Processed {total_reads} reads...", end='')

    print(f"\r    Processed {total_reads} reads. -- Done.")

    plot_mapped_percentages(counts, total_reads, read_length)
    plot_phred_distributions(phreds, read_length)
    bam.close()


def plot_mapped_percentages(
        counts: list[int], total_reads: int, k: int, fig_filename: "os.PathLike[str] | str" = 'map_perc.png'
) -> None:
    """Calculate and plot the percentage of mapped calls by read position.

    Args:
        counts (list[int]): A list of counts of mapped calls for each read position.
        total_reads (int): The total number of reads.
        k (int): The length of the reads.
        fig_filename (os.PathLike[str] | str): The filename to save the plot to.
    """
    print('Calculating mapped percentage frequencies... ', end='')
    freqs = [100 * x / total_reads for x in counts]
    print('Done.')

    print("Creating mapped percentage plot... ", end='')
    fig, ax = plt.subplots(figsize=(16, 9), dpi=300, tight_layout=True)
    ax.plot(range(1, k + 1), freqs)
    ax.set_xlabel('Read distance', fontsize='xx-large')
    ax.set_ylabel('Percentage of mapped calls', fontsize='xx-large')
    fig.suptitle(
        'Percentage of mapped calls as a function of the position from the start of the sequencer read',
        fontsize='xx-large'
    )
    fig.savefig(fig_filename)
    plt.close(fig)  # Close the figure to free memory
    print(f'\rSaved mapped percentage plot to "{fig_filename}"')


def plot_phred_distributions(
        phreds: dict[int, list[int]], k: int, fig_filename: "os.PathLike[str] | str" = 'phred2.png'
) -> None:
    """Calculate and plot the distribution of PHRED scores by read position.

    Args:
        phreds (dict[int, list[int]]): A dictionary with read positions as keys and lists of PHRED scores as values.
        k (int): The length of the reads.
        fig_filename (os.PathLike[str] | str): The path to save the plot to.
    """
    print(f'Calculating PHRED scores statistics for each position (0 - {k - 1}): ')
    print('  maxs', end='')
    maxs = [max(phreds[i]) for i in range(k)]
    print(', tops', end='')
    tops = [np.percentile(phreds[i], 95) for i in range(k)]
    print(', medians', end='')
    medians = [np.percentile(phreds[i], 50) for i in range(k)]
    print(', bottoms', end='... ')
    bottoms = [np.percentile(phreds[i], 5) for i in range(k)]
    print('Done.')

    print('Calculating differences for stackplot layers: ')
    print('  fig_medians', end='')
    fig_medians = [med - bot for med, bot in zip(medians, bottoms)]
    print(', fig_tops', end='')
    fig_tops = [top - med for top, med in zip(tops, medians)]
    print(', fig_maxs', end='... ')
    fig_maxs = [ma - top for ma, top in zip(maxs, tops)]
    print('Done.')

    print("Creating PHRED score distribution plot... ", end='')
    fig, ax = plt.subplots(figsize=(16, 9), dpi=300, tight_layout=True)
    ax.stackplot(range(1, k + 1), (bottoms, fig_medians, fig_tops, fig_maxs))
    ax.plot(range(1, k + 1), maxs, 'k-')  # Plot maximum scores as a black line
    ax.set_xlabel('Read distance', fontsize='xx-large')
    ax.set_ylabel('PHRED score', fontsize='xx-large')
    fig.suptitle(
        'Distribution of PHRED scores as a function of the position in the read',
        fontsize='xx-large'
    )
    fig.savefig(fig_filename)
    plt.close(fig)  # Close the figure to free memory
    print(f'\rSaved PHRED score distribution plot to "{fig_filename}"')


if __name__ == '__main__':
    # Set up command-line argument parsing
    a_parser = argparse.ArgumentParser(
        description="Calculates and plots "
                    "the percentage of mapped calls by read position "
                    "and "
                    "the distribution of PHRED scores by read position"
    )
    a_parser.add_argument(
        '-c', '--chromosome',
        default=20,
        help="Chromosome number to process (default: 20)"
    )
    a_parser.add_argument(
        '-f', '--fromurl',
        default="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/exome_alignment/",
        help="Base of the URL for downloading data "
             "(default: "
             "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/exome_alignment/"
             ")"
    )
    a_parser.add_argument(
        '-p', '--path',
        default='NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam',
        help="Path to save the file, filename must be the same as original "
             "(default: NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam)"
    )
    a_parser.add_argument(
        '-r', '--length',
        default=76,
        help="The length of the reads (default: 76)"
    )

    args = a_parser.parse_args()

    main(
        chromosome_number=args.chromosome,
        base_url=args.fromurl,
        bam_path=Path(args.path),
        read_length=args.length
    )
