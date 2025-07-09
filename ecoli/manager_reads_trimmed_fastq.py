#!/usr/bin/env python3
import argparse
import glob

import working_with_fastq as worfas


def main(pattern):
    """
    Performs FASTQ (gunzipped) per-file quality control and analysis
    of all files matching the given pattern.

    Args:
        pattern (str): The glob pattern as a string to match files against.

    Returns:
        None
    """
    msg = f"Processing all FASTQ (gunzipped) files matching pattern: {pattern} "
    print(msg)
    print("#" * len(msg[:-1]))
    for filename_fastq_gz in glob.glob(pattern):
        try:
            worfas.process_fastq_perform_analyses(filename_fastq_gz)
        except Exception as exn:
            print()
            print(f"! EXCEPTION occured with {filename_fastq_gz}: {exn}")
            print()


if __name__ == "__main__":
    pattern_default = "*.trim.fastq.gz"
    a_parser = argparse.ArgumentParser(description="Work with all FASTQ (gunzipped) files matching pattern")
    a_parser.add_argument(
        "--pattern", "-p",
        type=str,
        help=f"The pattern to match. Default: {pattern_default}",
        default=pattern_default,
    )
    args = a_parser.parse_args()
    main(args.pattern)
