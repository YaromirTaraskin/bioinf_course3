#!/usr/bin/env python3
""" (reads_1.fastq, reads_2.fastq |-> reference_genome.fasta) ==> aligned.sam """
import argparse
import os
import glob
import subprocess
from pathlib import Path


def main(dir_input, dir_output, genome_ref, suffix_fastq_trimmed, input_gzipped, path_bwa="bwa"):
    """ Main function that aligns paired-end FASTQ files to a reference genome using BWA.

    Args:
        dir_input (str | os.PathLike): Directory containing trimmed FASTQ files to align.
        dir_output (str | os.PathLike): Directory to store resulting SAM alignment files.
        genome_ref (str | os.PathLike): Path to reference genome FASTA file.
        suffix_fastq_trimmed (str): Suffix for trimmed FASTQ files.
        input_gzipped (bool | None): pass True to use gzipped input files or False value otherwise.
        path_bwa (str | os.PathLike): Path to BWA executable.
    """
    input_extension = ".fastq.gz" if input_gzipped else ".fastq"
    print(f'Working with files ending "{suffix_fastq_trimmed}{input_extension}" :')
    print(f'  chosen input extension is "{input_extension}" ')

    os.makedirs(dir_output, exist_ok=True)  # Ensure output directory exists

    print()
    print(f'Indexing reference genome "{genome_ref}" : ')
    subprocess.run(
        [str(arg) for arg in [path_bwa, "index", genome_ref]]
    )

    # Process all paired-end FASTQ files
    for path_fq1 in glob.glob(f"*_1{suffix_fastq_trimmed}{input_extension}", root_dir=dir_input):
        base_name = Path(path_fq1).stem.replace(f"_1{suffix_fastq_trimmed}", "")

        path_fq2 = Path(dir_input, f"{base_name}_2{suffix_fastq_trimmed}{input_extension}")
        path_sam = Path(dir_output, f"{base_name}.aligned.sam")

        print()
        print(f'Working with file "{base_name}" : ')

        try:
            # Run BWA alignment
            with open(path_sam, "w") as sam_file:
                subprocess.run(
                    [str(arg) for arg in [path_bwa, "mem", genome_ref, path_fq1, path_fq2]],
                    stdout=sam_file,
                    check=True
                )
        except subprocess.CalledProcessError as exn:
            print(f"Error occurred during alignment for {base_name}: {exn} ")
            if os.path.exists(path_sam):
                print(f"Removing incomplete alignment file: {path_sam} ")
                os.remove(path_sam)
        except Exception as exn:
            print(f"Unexpected error occurred during alignment for {base_name}: {exn} ")
        finally:
            print(f'Finished working with "{base_name}" .')


if __name__ == "__main__":
    a_parser = argparse.ArgumentParser(
        description='Align paired-end FASTQ files '
                    '(Pair format: "<name>_1<suffix>.fastq" and "<name>_2<suffix>.fastq") '
                    'to a reference genome (FASTA)'
                    ' into SAM file using BWA.'
    )
    a_parser.add_argument(
        "--dir-input", "-i",
        type=str,
        help="Directory containing trimmed FASTQ files to align. Default: current directory.",
        default="."
    )
    a_parser.add_argument(
        "--dir-output", "-o",
        type=str,
        help="Directory to store resulting SAM alignment files. Default: current directory.",
        default="."
    )
    a_parser.add_argument(
        "--genome-ref", "-r",
        type=str,
        help="Path to reference genome FASTA file. Default: ./ecoli_rel606.fasta",
        default="./ecoli_rel606.fasta"
    )
    a_parser.add_argument(
        "--suffix-fastq-trimmed", "-x",
        type=str,
        help="Suffix for trimmed FASTQ files. Default: .trim.sub",
        default=".trim.sub"
    )
    a_parser.add_argument(
        "--gunzipped", "-z",
        help='Flag to indicate that the input files are gunzipped (ends with ".gz").'
    )

    args = a_parser.parse_args()
    main(
        dir_input=args.dir_input,
        dir_output=args.dir_output,
        genome_ref=args.genome_ref,
        suffix_fastq_trimmed=args.suffix_fastq_trimmed,
        input_gzipped=args.gunzipped
    )
