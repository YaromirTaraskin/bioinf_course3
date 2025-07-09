#!/usr/bin/env python3

"""
A Python console program to analyze FASTQ files.
- Counts occurrences and percentages of each letter in sequences.
- Plots the number of 'N' calls per position.
- Counts occurrences and percentages of PHRED quality scores for positions >=26.
- Plots the distribution of PHRED scores per position (>=26, excluding 40).
"""
import os
import subprocess
from collections import defaultdict
import gzip
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO


def main():
    """Main function to process the FASTQ file and perform analyses."""
    # Check command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python3 Working_with_FASTQ.py <fastq.gz>")
        print("#")
        filename = "SRR003265.filt.fastq.gz"
        print(f'Using "{filename}" as example')
        wget_download_file_if_not_exists(
            "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz",
            filename
        )
    else:
        filename = sys.argv[1]

    try:
        # Open and parse the FASTQ file
        with gzip.open(filename, "rt", encoding="utf-8") as f:
            records = SeqIO.parse(f, "fastq")

            # Initialize data structures
            counts_letter = defaultdict(int)  # Letter counts
            counts_n_per_pos = defaultdict(int)  # 'N' counts per position
            quality_counts = defaultdict(int)  # Quality score counts
            quality_per_pos = defaultdict(list)  # Quality scores per position
            len_max = 0  # Maximum sequence length

            # Process each record in a single pass
            for i, rec in enumerate(records, start=1):
                print(f"\r  processing record {i} ... ", end='')
                sequence = rec.seq
                annotation_quality = rec.letter_annotations["phred_quality"]
                len_max = max(len_max, len(sequence))

                # Count letters in sequence
                for letter in sequence:
                    counts_letter[letter] += 1

                # Count 'N's per position
                for pos0, letter in enumerate(sequence):
                    if letter == "N":
                        counts_n_per_pos[pos0+1] += 1

                # Count quality scores and collect for distribution
                for pos0, qual in enumerate(annotation_quality):
                    if pos0 >= 25:  # Positions >= 26 (0-based index 25)
                        quality_counts[qual] += 1
                        if qual != 40:
                            quality_per_pos[pos0+1].append(qual)
            print("Done.")
            # Analyze and output results
            analyze_letters(counts_letter)
            plot_n_calls(counts_n_per_pos, "n_calls.png")
            analyze_quality_scores(quality_counts)
            plot_quality_distributions(quality_per_pos, "phred.png")

    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        sys.exit(1)
    except Exception as exn:
        print(f"Error: An unexpected error occurred: {exn}")
        sys.exit(1)


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


def analyze_letters(letter_counter: "dict[str, int]") -> None:
    """Calculate and print letter counts and percentages.

    Args:
        letter_counter (dict[str, int]): A dictionary with letters as keys and their counts as values.
    """
    total = sum(letter_counter.values())
    if total == 0:
        print("No sequences found.")
        return

    print("Letter counts:")
    for letter, l_count in sorted(letter_counter.items()):
        l_percentage = 100 * l_count / total
        print(f"  {letter}: {l_percentage:.2f}% {l_count}")


def plot_n_calls(counter_ns: "dict[int, int]", fig_filename: "os.PathLike[str] | str") -> None:
    """Plot the number of 'N' calls per position and save to file.

    Args:
        counter_ns: A dictionary with read positions as keys and counts of 'N' calls as values.
        fig_filename: The filename to save the plot to.
    """
    if not counter_ns:
        print("No 'N' calls found to plot.")
        return
    print("Plotting number of 'N' calls per position: ", end="")
    seq_len = max(counter_ns.keys())
    positions = range(1, seq_len + 1)
    n_values = [counter_ns[x] for x in positions]

    fig, ax = plt.subplots(figsize=(16, 9), dpi=300)
    ax.plot(positions, n_values)
    ax.set_xlim(1, seq_len)
    ax.set_xlabel("Read distance", fontsize="xx-large")
    ax.set_ylabel("Number of N Calls", fontsize="xx-large")
    fig.suptitle(
        "Number of N calls as a function of the distance from the start of the sequencer read",
        fontsize="xx-large",
    )
    fig.tight_layout()
    fig.savefig(fig_filename)
    plt.close(fig)
    print(f"\rSaved 'N' calls plot to '{fig_filename}'")


def analyze_quality_scores(counter_quality: "dict[int, int]") -> None:
    """Calculate and print quality score counts and percentages.

    Args:
        counter_quality: A dictionary with quality scores as keys and counts as values.
    """
    total_quality = sum(counter_quality.values())
    if total_quality == 0:
        print("\nNo quality scores found for positions >=26.")
        return

    print("\nQuality score counts (positions >=26):")
    for qual, q_count in sorted(counter_quality.items()):
        q_percentage = 100 * q_count / total_quality
        print(f"  {qual}: {q_percentage:.2f}% {q_count}")


def plot_quality_distributions(
        quality_per_position: "dict[int, list[int]]", fig_filename: "os.PathLike[str] | str"
) -> None:
    """Plot the distribution of quality scores per position and save to file.

    Args:
        quality_per_position: A dictionary with read positions as keys and lists of quality scores as values.
        fig_filename: The path to save the plot to.
    """
    if not quality_per_position:
        print("No quality scores collected for plotting (positions >=26, excluding 40).")
        return

    print("Preparing data for boxplot of quality scores per position: ", end="")

    positions = sorted(quality_per_position.keys())
    value_per_pos = [quality_per_position[pos] for pos in positions]

    print(f"\rPlotting quality score distributions... ", end="")

    fig, ax = plt.subplots(figsize=(16, 9), dpi=300)
    sns.boxplot(data=value_per_pos, ax=ax)
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels([str(pos) for pos in positions])
    ax.set_xlabel("Read distance", fontsize="xx-large")
    ax.set_ylabel("PHRED score", fontsize="xx-large")
    fig.suptitle(
        "Distribution of PHRED scores as a function of read distance (positions >=26, excluding 40)",
        fontsize="xx-large",
    )
    fig.tight_layout()
    fig.savefig(fig_filename)
    plt.close(fig)

    print(f"\rSaved PHRED score distribution plot to '{fig_filename}'")


if __name__ == "__main__":
    main()
