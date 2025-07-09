#!/usr/bin/env python3

"""Analyze a BED file using HTSeq."""

import argparse
from collections import defaultdict
import re
from io import StringIO
import HTSeq


def main(bed_file):
    """Analyze a BED file using HTSeq.

    Args:
        bed_file (os.PathLike | str): Path to the BED file
    """
    # Step 1: Read the BED file and filter out lines with fewer than 3 fields
    with open(bed_file, 'r') as b_file:
        bed_lines_valid = [line for line in b_file if len(line.strip().split('\t')) >= 3]
    bed_content = StringIO(''.join(bed_lines_valid))
    bed_lct = HTSeq.BED_Reader(bed_content)

    # Step 2: Initialize data structures
    feature_types = defaultdict(int)  # Counts of feature types (e.g., ENSE, CCDS)
    exon_start = None                 # Min start position for CCDS records
    exon_end = None                   # Max end position for CCDS records
    sizes = []                        # List of CCDS exon lengths
    length_total = 0                  # Sum of CCDS exon lengths for mean
    count = 0                         # Number of CCDS records
    rec_last = None                   # Last processed record

    # Step 3: Process each record in a single pass
    for record in bed_lct:
        rec_last = record  # Update last record

        # Extract feature type (uppercase letters at start of name)
        match = re.search('([A-Z]+)', record.name)
        if match:
            feature_type = match.group(0)
            feature_types[feature_type] += 1
        else:
            print(f"Warning: No feature type found in record name: {record.name}")

        # Process CCDS records
        if record.name.startswith('CCDS'):
            interval = record.iv
            # Update overall start and end
            if exon_start is None:
                exon_start = interval.start
            else:
                exon_start = min(exon_start, interval.start)
            if exon_end is None:
                exon_end = interval.end
            else:
                exon_end = max(exon_end, interval.end)
            # Compute and store size
            size = interval.length
            sizes.append(size)
            length_total += size
            count += 1

    # Step 4: Output results
    # Feature type counts
    print("Feature Types:")
    for ft_name, ft_count in feature_types.items():
        print(f"  {ft_name}:\t{ft_count}")

    # Last record details
    if rec_last:
        print("\nLast Record Details:")
        print(f"  Record: {rec_last}")
        print(f"  Name: {rec_last.name}")
        print(f"  Type: {type(rec_last)}")
        interval = rec_last.iv
        print(f"  Interval: {interval}")
        print(f"  Type of interval: {type(interval)}")
        print(f"  Chromosome: {interval.chrom}")
        print(f"  Start: {interval.start}")
        print(f"  End: {interval.end}")
        print(f"  Strand: {interval.strand}")
        print(f"  Length: {interval.length}")
        print(f"  Start_d: {interval.start_d}")
        print(f"  Start as position: {interval.start_as_pos}")

    # CCDS exon statistics
    print("\nExon Statistics -- for records starting with 'CCDS':")
    if count > 0:
        size_min = min(sizes)
        size_max = max(sizes)
        size_mean = length_total / count
        print(f"  Number of exons: {count}")
        print(f"  Start position: {exon_start}")
        print(f"  End position: {exon_end}")
        print(f"  Smallest exon size: {size_min}")
        print(f"  Largest exon size: {size_max}")
        print(f"  Mean exon size: {size_mean:.1f}")
    else:
        print("  No exons with names starting with 'CCDS' were found.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze a BED file using HTSeq.")
    parser.add_argument(
        '-f', '--bedfile',
        default='LCT.bed',
        help='Path to the BED file to analyze, default: "LCT.bed"'
    )
    args = parser.parse_args()
    main(args.bedfile)
