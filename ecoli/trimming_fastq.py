#!/usr/bin/env python3
"""
Trimmomatic main functionality Python emulation (based on Biopython)

Engineered by Yaromir Taraskin
"""


import argparse
import glob
import gzip

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_adapters(fasta_file):
    """Load adapter sequences and their names from a FASTA file.

    Args:
        fasta_file (str | os.PathLike): Path to the FASTA file containing adapter sequences.

    Returns:
        list[tuple[str, str]]: A list of tuples containing adapter names and sequences.
    """
    adapters_named = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        adapters_named.append((record.id, str(record.seq)))

    return adapters_named


def find_seed_positions(read_seq, seed, mismatches_max):
    """Find positions where the seed matches the read scanned from 3' end.

    Args:
        read_seq (str): The read sequence.
        seed (str): The seed sequence.
        mismatches_max (int): The maximum number of mismatches allowed in the seed region.

    Returns:
        list[int]: A list of positions where the seed matches the read.
    """
    len_seed = len(seed)
    len_read = len(read_seq)
    positions = []

    for i in range(len_read - len_seed, -1, -1):
        mismatches = 0
        for j in range(len_seed):
            if read_seq[i + j] != seed[j]:
                mismatches += 1
                if mismatches > mismatches_max:
                    break
        if mismatches <= mismatches_max:
            positions.append(i)

    return positions


def calculate_score(read_seq, adapter, pos_start):
    """Calculate alignment score (number of matching bases).

    Args:
        read_seq (str): The read sequence.
        adapter (str): The adapter sequence.
        pos_start (int): The starting position for the comparison with the adapter.

    Returns:
        int: The number of matching bases.
    """
    score = 0
    adapter_len = len(adapter)
    read_len = len(read_seq)

    for i in range(adapter_len):
        if pos_start + i >= read_len:
            break
        if read_seq[pos_start + i] == adapter[i]:
            score += 1

    return score


def find_adapter_matches(read, adapters_named, seed_mismatches, *, seed_size=16):
    """Find all adapter matches in a read.

    Args:
        read (SeqRecord): The read to search for adapter matches.
        adapters_named (list[tuple[str, str]]):
            A list of tuples containing adapter names and sequences.
        seed_mismatches (int): The maximum number of mismatches allowed in the seed region.
        seed_size (int): The length of the seed region.

    Returns:
        list[tuple[str, str, int, int]]:
            A list of tuples containing adapter names, sequences, positions, and scores.
    """
    read_seq = str(read.seq)
    matches = []

    for a_name, adapter in adapters_named:
        if len(adapter) < seed_size:
            continue  # Skip short adapters
        seed = adapter[:seed_size]
        positions = find_seed_positions(read_seq, seed, seed_mismatches)
        for pos in positions:
            score = calculate_score(read_seq, adapter, pos)
            matches.append((a_name, adapter, pos, score))

    return matches


def process_individual(read, matches, threshold_clip_simple):
    """Trim a read based on the best individual adapter match.

    Args:
        read (SeqRecord): The read to process.
        matches (list[tuple[str, str, int, int]]):
            A list of tuples containing adapter names, sequences, positions, and scores.
        threshold_clip_simple (int): The threshold for simple clipping.

    Returns:
        tuple(SeqRecord, str | None, int):
            A tuple containing the trimmed read, adapter name, and the alignment score.
    """
    score_best = -1
    pos_best = None
    name_best = None

    for name_cur, _adapter, pos_cur, score_cur in matches:
        if score_cur > score_best and score_cur >= threshold_clip_simple:
            name_best, pos_best, score_best = name_cur, pos_cur, score_cur

    if score_best >= threshold_clip_simple:
        return read[:pos_best], name_best, score_best

    return read, None, 0


def sliding_window_trim(read, window_size, quality_avg_min):
    """
    Apply SLIDING WINDOW quality trimming to a read.

    Args:
        read (SeqRecord): Input read with sequence and quality scores.
        window_size (int): Size of the sliding window (bases).
        quality_avg_min (int): Minimum average quality required in the window.

    Returns:
        SeqRecord: Trimmed read.
    """
    seq = str(read.seq)
    quals = read.letter_annotations.get('phred_quality', [])

    if not quals or len(seq) != len(quals):
        return read  # Skip if no quality data
    if len(seq) < window_size:
        return read  # Too short to process

    pos_trim = None

    for pos_cur in range(len(seq) - window_size + 1):
        window = quals[pos_cur: pos_cur + window_size]
        quality_avg_cur = sum(window) / len(window)
        if quality_avg_cur < quality_avg_min:
            pos_trim = pos_cur
            break  # Trim at the first failing window

    if pos_trim is None:
        return read  # No trimming needed

    # Trim to the start of the failing window
    return read[:pos_trim]


def trim_pair(
        read_1, read_2, adapters,
        window_size, threshold_qual,
        len_min,
        mismatches_seed, threshold_clip_palindrome, threshold_clip_simple
):
    """
    Trim adapter sequences from a pair of reads.

    Args:
        read_1 (SeqRecord): Forward read.
        read_2 (SeqRecord): Reverse read.
        adapters (list): List of (name, sequence) adapters.
        mismatches_seed (int): Allowed seed mismatches.
        threshold_clip_palindrome (int): Score threshold for palindrome clipping.
        threshold_clip_simple (int): Score threshold for simple clipping.
        len_min (int): Minimum length to retain reads.
        window_size (int): The size of the sliding window for quality trimming.
        threshold_qual (int): The minimum mean quality score for a window.

    Returns:
        tuple: (paired_reads, unpaired_reads, discarded_reads)
    """
    # Find all adapter matches for both reads
    read_1_matches = find_adapter_matches(read_1, adapters, mismatches_seed)
    read_2_matches = find_adapter_matches(read_2, adapters, mismatches_seed)

    palindrome_best = None
    palindrome_best_score = -1
    # Check for palindrome pairs -- reverse-complement adapters
    for r1_name, r1_adapter, r1_pos, r1_score in read_1_matches:
        for r2_name, r2_adapter, r2_pos, r2_score in read_2_matches:
            if (
                    r2_name == f'{r1_name}_rc'
                    or
                    r1_name == f'{r2_name}_rc'
            ):
                combined_score = r1_score + r2_score
                if (
                        combined_score > palindrome_best_score
                        and
                        combined_score >= threshold_clip_palindrome
                ):
                    palindrome_best = (r1_pos, r2_pos, combined_score)
                    palindrome_best_score = combined_score

    # Trim based on palindrome or individual matches
    if palindrome_best:
        r1_pos_trim, r2_pos_trim, _ = palindrome_best
        trimmed_1 = read_1[:r1_pos_trim]
        trimmed_2 = read_2[:r2_pos_trim]
    else:
        trimmed_1 = process_individual(read_1, read_1_matches, threshold_clip_simple)[0]
        trimmed_2 = process_individual(read_2, read_2_matches, threshold_clip_simple)[0]

    # SLIDINGWINDOW quality trimming
    trimmed_1 = sliding_window_trim(trimmed_1, window_size, threshold_qual)
    trimmed_2 = sliding_window_trim(trimmed_2, window_size, threshold_qual)

    # MINLEN filter
    paired_1 = trimmed_1 if len(trimmed_1) >= len_min else None
    paired_2 = trimmed_2 if len(trimmed_2) >= len_min else None

    unpaired_1 = trimmed_1 if (paired_1 and not paired_2) else None
    unpaired_2 = trimmed_2 if (paired_2 and not paired_1) else None

    return (paired_1, paired_2), (unpaired_1, unpaired_2)


def process_file_pair(
        infile_base, adapter_file,
        window_size, threshold_qual,
        len_min,
        mismatches_seed, clip_palindrome, clip_simple
):
    """
    Process paired-end FASTQ files, perform trimming, and write trimmed sequences to output files.

    This function reads paired-end FASTQ files,
    applies adapter and quality trimming, and
    writes the trimmed sequences to new FASTQ files.

    It handles paired reads and writes paired trimming results to
    output files with ".trim" suffix.
    It also handles reads got unpaired after trimming and writes them to
    separate output files with "un.trim" suffix.

    Args:
        infile_base (str | os.PathLike):
            Base name of the input FASTQ files (without "_1.fastq.gz" or "_2.fastq.gz" suffix).
        adapter_file (str | os.PathLike): Path to the adapters_named file.
        window_size (int): The size of the sliding window for quality trimming.
        threshold_qual (int): The minimum mean quality score for a window.
        len_min (int): The minimum length of the sequence to be considered.
        mismatches_seed (int): The number of mismatches allowed in the seed region.
        clip_palindrome (int): Score threshold for palindrome clipping.
        clip_simple (int): Score threshold for simple clipping.

    Returns:
        None
    """
    print(f'Processing "{infile_base}" : ')
    adapters_named = load_adapters(adapter_file)

    # Initialize counters for statistics
    cnt_paired = 0
    cnt_unpaired_fwd = 0
    cnt_unpaired_rev = 0

    # Open input and output FASTQ files
    with gzip.open(f"{infile_base}_1.fastq.gz", mode="rt") as infile_fwd, \
            gzip.open(f"{infile_base}_2.fastq.gz", mode="rt") as infile_rev, \
            gzip.open(f"{infile_base}_1.trim.fastq.gz", mode="wt") as hdl_fwd_paired, \
            gzip.open(f"{infile_base}_1un.trim.fastq.gz", mode="wt") as hdl_fwd_unpaired, \
            gzip.open(f"{infile_base}_2.trim.fastq.gz", mode="wt") as hdl_rev_paired, \
            gzip.open(f"{infile_base}_2un.trim.fastq.gz", mode="wt") as hdl_rev_unpaired:

        # Iterate through each pair of reads
        for idx, (fwd, rev) in enumerate(
                zip(
                    SeqIO.parse(infile_fwd, format="fastq"),
                    SeqIO.parse(infile_rev, format="fastq")
                ),
                start=1
        ):
            print(f"\r  processing records pair {idx} ... ", end='')

            # Perform trimming on the read pair
            (paired_1, paired_2), (unpaired_1, unpaired_2) = trim_pair(
                read_1=fwd, read_2=rev,
                adapters=adapters_named,
                window_size=window_size, threshold_qual=threshold_qual,
                len_min=len_min,
                mismatches_seed=mismatches_seed,
                threshold_clip_palindrome=clip_palindrome,
                threshold_clip_simple=clip_simple
            )

            # Write trimmed output to respective FASTQ files
            if paired_1 and paired_2:
                SeqIO.write(paired_1, hdl_fwd_paired, format="fastq")
                SeqIO.write(paired_2, hdl_rev_paired, format="fastq")
                cnt_paired += 1
            else:
                if unpaired_1:
                    SeqIO.write(unpaired_1, hdl_fwd_unpaired, format="fastq")
                    cnt_unpaired_fwd += 1

                if unpaired_2:
                    SeqIO.write(unpaired_2, hdl_rev_unpaired, format="fastq")
                    cnt_unpaired_rev += 1
    print("Done.")
    print(f"Stats:\n"
          f"  Kept pairs: {cnt_paired}\n"
          f"  Kept unpaired forward: {cnt_unpaired_fwd}\n"
          f"  Kept unpaired reverse: {cnt_unpaired_rev}\n")
    print()


def main(
        adapter_file,
        window_size, threshold_qual,
        len_min,
        mismatches_seed, clip_palindrome, clip_simple
):
    """
    Main function to trim a set of FASTQ files in the current working directory.

    This function takes the following parameters:
        adapter_file (str | os.PathLike): The path to a file containing the adapter sequences.
        window_size (int): The size of the sliding window for quality trimming.
        threshold_qual (int): The minimum mean quality score for a window.
        len_min (int): The minimum length of the sequence to be considered.
        mismatches_seed (int): The number of mismatches allowed in the seed region.
        clip_palindrome_perc (int): The percentage of the read length to be clipped
            from the palindrome region.
        clip_simple_perc (int): The percentage of the read length to be clipped
            from the simple region.

    The function will trim the FASTQ files
    in the current working directory using
    the specified parameters and
    write the trimmed sequences to new FASTQ files
    with the same base name as the original file but
    with "_1.trim.fastq.gz" and "_2.trim.fastq.gz" appended to the end
    for the trimmed pairs.

    The function will also
    write the trimmed sequences that became unpaired after trimming
    to new FASTQ files with the same base name as the original file but
    with "_1un.trim.fastq.gz" and "_2un.trim.fastq.gz" appended to the end.

    The function does not return anything.
    """
    for infile in glob.glob("*_1.fastq.gz"):
        base = infile.replace("_1.fastq.gz", "")
        process_file_pair(
            infile_base=base, adapter_file=adapter_file,
            window_size=window_size, threshold_qual=threshold_qual,
            len_min=len_min,
            mismatches_seed=mismatches_seed, clip_palindrome=clip_palindrome, clip_simple=clip_simple
        )


if __name__ == "__main__":
    a_parser = argparse.ArgumentParser(
        description='Trim paired-end FASTQ files using Trimmomatic analog written in Python.'
    )
    a_parser.add_argument(
        "--adapter-file", "-a",
        help='File containing adapter sequences. Default: "NexteraPE-PE.fa"',
        default="NexteraPE-PE.fa"
    )
    a_parser.add_argument(
        "--window-size", "-w",
        help='Size of the sliding window for quality trimming. Default: 4',
        default=4
    )
    a_parser.add_argument(
        "--threshold-qual", "-q",
        help='Minimum mean quality score for a sliding window. Default: 20',
        default=20
    )
    a_parser.add_argument(
        "--len-min", "-l",
        help='Minimum length of the sequence to be accepted after trimming. Default: 25',
        default=25
    )
    a_parser.add_argument(
        "--mismatches-seed", "-m",
        help='Maximum number of mismatches allowed in the seed region. Default: 2',
        default=2
    )
    a_parser.add_argument(
        "--clip-palindrome", "-p",
        help='Threshold for palindrome clipping: '
             'the sum of the alignment scores of both reads for adapter-adapter hybrids '
             '(dimers) where both reads in a pair contain adapters. Default: 40',
        default=40
    )
    a_parser.add_argument(
        "--clip-simple", "-s",
        help='Threshold for simple clipping: individual alignment score for non-dimer adapters. '
             'Default: 15',
        default=15
    )

    args = a_parser.parse_args()
    main(
        adapter_file=args.adapter_file,
        window_size=args.window_size,
        threshold_qual=args.threshold_qual,
        len_min=args.len_min,
        mismatches_seed=args.mismatches_seed,
        clip_palindrome=args.clip_palindrome,
        clip_simple=args.clip_simple
    )
