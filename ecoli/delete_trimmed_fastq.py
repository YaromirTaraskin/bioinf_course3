#!/usr/bin/env python3
import argparse
import glob
import os


def main(pattern):
    """
    Deletes all files matching the given pattern.

    Args:
        pattern (str): The glob pattern as a string to match files against.

    Returns:
        None
    """
    msg = f"Deleting all files matching pattern: {pattern} "
    print(msg)
    print("=" * len(msg[:-1]))
    print()
    for filename in glob.glob(pattern):
        print(f'Deleting "{filename}" ... ', end="")
        try:
            os.remove(filename)
            print("Done. ")
        except OSError as os_err:
            print(f" ! OS ERROR : {os_err}")
        except Exception as exn:
            print(f" ! UNEXPETED EXCEPTION : {exn}")
        print()


if __name__ == "__main__":
    pattern_default = "*.trim.fastq.gz"
    a_parser = argparse.ArgumentParser(description="Delete all files matching a pattern!")
    a_parser.add_argument(
        "--pattern", "-p",
        type=str,
        help=f"The pattern to match. Default: {pattern_default}",
        default=pattern_default,
    )
    args = a_parser.parse_args()
    main(args.pattern)
