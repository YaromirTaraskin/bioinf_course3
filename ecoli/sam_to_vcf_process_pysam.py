#!/usr/bin/env python3
"""
aligned.sam ==> .bam ==> sorted.bam ==> indexed sorted.bam ==>
==> genotype_likelihoods_raw.bcf ==> variants.vcf
"""
import argparse
import glob
import os
import sys
from pathlib import Path

import pysam
import pysam.bcftools


def process_samtools(file_sam, file_bam, file_bam_sorted):
    """
    Convert SAM to BAM, sort, and index using SAMtools. (Part 1)

    Args:
        file_sam (str | os.PathLike): Input SAM file path.
        file_bam (str | os.PathLike): Output BAM file path.
        file_bam_sorted (str | os.PathLike): Sorted BAM file path.
    """
    # Convert file paths to strings for using with pysam commands
    file_bam = str(file_bam)
    file_bam_sorted = str(file_bam_sorted)

    print(f'Converting SAM to BAM : "{file_sam}" ==> "{file_bam}" ... ', end='', flush=True)
    with pysam.AlignmentFile(file_sam, "r") as infile_sam:
        with pysam.AlignmentFile(file_bam, "wb", header=infile_sam.header) as outfile_bam:
            for read in infile_sam:
                outfile_bam.write(read)
    print("Done. ")

    print(f'Sorting BAM : "{file_bam}" ==> "{file_bam_sorted}" ... ', end="", flush=True)
    pysam.sort("-o", file_bam_sorted, file_bam)
    print("Done. ")

    print(f'Indexing sorted BAM : "{file_bam_sorted}" ... ', end="", flush=True)
    pysam.index(file_bam_sorted)
    print("Done. ")


def process_bcftools(genome_ref_fasta, bam_sorted, bcf_raw, vcf_variants):
    """
    Generate variant calls using BCFtools. (Part 2)

    Args:
        genome_ref_fasta (str | os.PathLike): Reference genome FASTA file path.
        bam_sorted (str | os.PathLike): Sorted BAM file path.
        bcf_raw (str | os.PathLike): Raw BCF output file path.
        vcf_variants (str | os.PathLike): Final variants VCF file path.
    """
    # Convert file paths to strings for using with pysam.bcftools commands
    genome_ref_fasta, bam_sorted, bcf_raw, vcf_variants = [
        str(arg) for arg in (genome_ref_fasta, bam_sorted, bcf_raw, vcf_variants)
    ]

    print(
        f'Generating raw BCF : ("{bam_sorted}", "{genome_ref_fasta}") ==> "{bcf_raw}" ... ',
        end="", flush=True
    )
    pysam.bcftools.mpileup(
        "-O", "b",  # Output in compressed BCF format
        "-o", bcf_raw,  # output to file, rather than stdout
        "--fasta-ref", genome_ref_fasta,  # faidx-indexed reference file in the FASTA format
        bam_sorted,  # from
        catch_stdout=False  # don't capture stdout, required for writing to file
    )
    print("Done. ")

    print(f'Generating variants VCF : "{bcf_raw}" ==> "{vcf_variants}" ... ', end='', flush=True)
    pysam.bcftools.call(
        "--ploidy", "1",  # Haploid organism
        "-m",  # Use multi-allelic caller
        "-v",  # Output variants only
        "-o", vcf_variants,  # write single stream output to file rather than stdout
        bcf_raw,  # from
        catch_stdout=False  # don't capture stdout, required for writing to file
    )
    print("Done. ")


def work_with_filenames_single_pipe(
        path_genome_ref,
        path_sam_input,
        path_samtools_bam,
        path_samtools_bam_sorted,
        path_bcftools_bcf_raw,
        path_vcf_output
):
    """
    Perform a single pipe
    (meaning only one SAM input file and only one VCF output file)
    with samtools and bcftools.

    Args:
        path_genome_ref (str | os.PathLike): Reference genome FASTA file path.
        path_sam_input (str | os.PathLike): Input SAM file path.
        path_samtools_bam (str | os.PathLike):
            BAM middle output file path (will be produced from SAM).
        path_samtools_bam_sorted (str | os.PathLike):
            Sorted BAM middle output file path (will be produced from BAM).
        path_bcftools_bcf_raw (str | os.PathLike):
            genotype likelihoods for alignments raw BCF middle output file path
            (will be produced from sorted BAM).
        path_vcf_output (str | os.PathLike):
            variants VCF output file path (will be produced from raw BCF).
    """
    print("samtools: ")
    print("---------")
    process_samtools(
        file_sam=path_sam_input,
        file_bam=path_samtools_bam,
        file_bam_sorted=path_samtools_bam_sorted
    )
    print()
    print("bcftools: ")
    print("---------")
    process_bcftools(
        genome_ref_fasta=path_genome_ref,
        bam_sorted=path_samtools_bam_sorted,
        bcf_raw=path_bcftools_bcf_raw,
        vcf_variants=path_vcf_output
    )
    print()


def main(
        path_genome_ref,
        suffix_sam_input,
        suffix_samtools_bam,
        suffix_samtools_bam_sorted,
        suffix_bcftools_bcf_raw,
        suffix_vcf_output,
        dir_input,
        dir_medium,
        dir_output
):
    """
    This function recognizes all SAM files with specified suffix
    in the specified input directory.

    For each input SAM file, this function
    will call `work_with_filenames_single_pipe`
    to process it according to the pipeline.

    If any exception occurs during processing,
    this function will print an error message
    to stderr and continue to the next file.

    Args:
        path_genome_ref (str | os.PathLike):
            the path to the *reference genome* FASTA file
        suffix_sam_input (str): the suffix for INPUT SAM files
        suffix_samtools_bam (str):
            the suffix for SAMtools BAM (middle output files)
        suffix_samtools_bam_sorted (str):
            the additional suffix for SAMtools sorted BAM (middle output files)
        suffix_bcftools_bcf_raw (str):
            the suffix for BCFtools raw BCF (middle output files)
        suffix_vcf_output (str): the suffix for VCF OUTPUT files
        dir_input (str | os.PathLike):
            the directory containing the input SAM files
        dir_medium (str | os.PathLike): the directory for middle output files
        dir_output (str | os.PathLike): the directory for VCF output files
    """
    os.makedirs(dir_medium, exist_ok=True)
    os.makedirs(dir_output, exist_ok=True)
    print(f'Using reference genome FASTA: "{path_genome_ref}" ')
    print()

    for basename_input in glob.glob(f"*{suffix_sam_input}.sam", root_dir=dir_input):
        path_sam_input = Path(dir_input, basename_input)
        base = Path(path_sam_input).stem
        base = base.removesuffix(suffix_sam_input)

        msg = f"Working with {base} : "
        print()
        print(msg)
        print("=" * len(msg[:-1]))
        print()

        path_samtools_bam = Path(
            dir_medium, f"{base}{suffix_samtools_bam}.bam"
        )
        path_samtools_bam_sorted = Path(
            dir_medium, f"{base}{suffix_samtools_bam_sorted}{suffix_samtools_bam}.bam"
        )
        path_bcftools_bcf_raw = Path(
            dir_medium, f"{base}{suffix_bcftools_bcf_raw}.bcf"
        )
        path_vcf_output = Path(
            dir_output, f"{base}{suffix_vcf_output}.vcf"
        )

        try:
            work_with_filenames_single_pipe(
                path_genome_ref=path_genome_ref,
                path_sam_input=path_sam_input,
                path_samtools_bam=path_samtools_bam,
                path_samtools_bam_sorted=path_samtools_bam_sorted,
                path_bcftools_bcf_raw=path_bcftools_bcf_raw,
                path_vcf_output=path_vcf_output
            )
        except Exception as exn:
            print("", file=sys.stderr)
            print(f"! EXCEPTION occurred with {base}: {exn}", file=sys.stderr)
            print("", file=sys.stderr)


if __name__ == "__main__":
    genome_default = "ecoli_rel606.fasta"
    a_parser = argparse.ArgumentParser(
        description="Convert SAM to BAM, sort, and index using SAMtools, "
                    "then generate variant calls using BCFtools."
    )
    a_parser.add_argument(
        "--genome-ref", "-g",
        type=str,
        help='Reference genome FASTA file path. '
             'Default: "ecoli_rel606.fasta" from --dir-input',
        default=None
    )
    a_parser.add_argument(
        "--suffix-sam-input", "-i",
        type=str,
        help="Suffix for INPUT alignment SAM file. Default: .aligned",
        default=".aligned"
    )
    a_parser.add_argument(
        "--suffix-bam", "-b",
        help='Suffix for SAMtools output BAM file. '
             'Default: same as --suffix-sam-input',
        default=None
    )
    a_parser.add_argument(
        "--suffix-bam-sorted", "-s",
        type=str,
        help="Additional suffix for sorted SAMtools output BAM file. "
             "Default: .sorted",
        default=".sorted"
    )
    a_parser.add_argument(
        "--suffix-bcf-raw", "-r",
        type=str,
        help="Suffix for BCFtools raw output BCF file. Default: _raw",
        default="_raw"
    )
    a_parser.add_argument(
        "--suffix-vcf-output", "-o",
        type=str,
        help="Suffix for OUTPUT variants VCF file. Default: _variants",
        default="_variants"
    )
    a_parser.add_argument(
        "--dir-input", "-d",
        type=str,
        help="Directory containing input SAM files. Default: current directory",
        default="."
    )
    a_parser.add_argument(
        "--dir-medium", "-m",
        help="Directory to store intermediate files. "
             "Default: same as --dir-output",
        default=None
    )
    a_parser.add_argument(
        "--dir-output", "-D",
        help="Directory to store output files. Default: same as --dir-input",
        default=None
    )

    args = a_parser.parse_args()
    if args.genome_ref is None:
        args.genome_ref = Path(args.dir_input, genome_default)

    main(
        path_genome_ref=args.genome_ref,
        dir_input=args.dir_input,
        dir_medium=args.dir_medium or args.dir_output or args.dir_input,
        dir_output=args.dir_output or args.dir_input,
        suffix_sam_input=args.suffix_sam_input,
        suffix_samtools_bam=args.suffix_bam or args.suffix_sam_input,
        suffix_samtools_bam_sorted=args.suffix_bam_sorted,
        suffix_bcftools_bcf_raw=args.suffix_bcf_raw,
        suffix_vcf_output=args.suffix_vcf_output
    )
