#!/usr/bin/env python3
"""
Filtering VCF files identified by filename suffix based on several quality metrics.

Constructed by Yaromir Taraskin, based on vcfutils.pl by lh3.
"""
import argparse
import glob
import os
import sys
from pathlib import Path

from cyvcf2 import VCF


def var_filter_outing(
    staging_item,
    is_print,
    output_main=sys.stdout,
    output_filtered=sys.stderr
):
    """(Former varFilter_aux) Prints a variant based on its filter status.

    Args:
        staging_item: List containing [score, flt_tag, indel_span, variant].
        is_print: Boolean indicating whether to print filtered variants to stderr.
        output_main: File object for main output.
        output_filtered: File object for filtered output.
    """
    score, flt_tag, indel_span, variant = staging_item
    if flt_tag == 0:
        print(str(variant), end='', file=output_main)  # Avoid extra newline
    elif is_print:
        tag = "UQdDaGgPMS"[flt_tag]
        print(f"{tag}\t{str(variant)}", file=output_filtered)


def parse_args(extended=False):
    """Parse command-line arguments.

    Args:
        extended (bool):  use the modified parameters list.
    """
    a_parser = argparse.ArgumentParser(description='Filter variants in VCF file.')
    a_parser.add_argument(
        '--quality-rms-min', '-Q',
        type=int,
        default=10, help='minimum RMS mapping quality for SNPs. Default: 10'
    )
    a_parser.add_argument(
        '--depth-read-min', '-d',
        type=int,
        default=2, help='minimum read depth. Default: 2'
    )
    a_parser.add_argument(
        '--depth-read-max', '-D',
        type=int,
        default=10000000, help='maximum read depth. Default: 10000000'
    )
    a_parser.add_argument(
        '--alternatives-min', '-a',
        type=int,
        default=2, help='minimum number of alternate bases. Default: 2'
    )
    a_parser.add_argument(
        '--wipe-int-around', '-w',
        type=int,
        default=3, help='SNP within INT bp around a gap to be filtered. Default: 3'
    )
    a_parser.add_argument(
        '--window-size', '-W',
        type=int,
        default=10, help='window size for filtering adjacent gaps. Default: 10'
    )
    a_parser.add_argument(
        '--p-min-strand', '-1',
        type=float,
        default=1e-4, help='min P-value for strand bias (given PV4). Default: 1e-4'
    )
    a_parser.add_argument(
        '--p-min-baseq', '-2',
        type=float,
        default=1e-100, help='min P-value for baseQ bias. Default: 1e-100'
    )
    a_parser.add_argument(
        '--p-min-mapq', '-3',
        type=float,
        default=0, help='min P-value for mapQ bias. Default: 0'
    )
    a_parser.add_argument(
        '--p-min-distance-end', '-4',
        type=float,
        default=1e-4, help='min P-value for end distance bias. Default: 1e-4'
    )
    a_parser.add_argument(
        '--p-min-hwe', '-e',
        type=float,
        default=1e-4, help='min P-value for HWE (plus F<0). Default: 1e-4'
    )
    a_parser.add_argument(
        '--mxgq-min', '-G',
        type=int,
        default=0, help='min MXGQ. Default: 0'
    )
    a_parser.add_argument(
        '--mxsp-max', '-S',
        type=int,
        default=1000, help='max MXSP. Default: 1000'
    )
    a_parser.add_argument(
        '--print-filtered', '-p',
        action='store_true',
        help='if set, print filtered variants'
    )
    if extended:
        a_parser.add_argument(
            "--suffix-input", "-i",
            type=str,
            help="Suffix for [identifying] input VCF files. Default: _variants",
            default="_variants"
        )
        a_parser.add_argument(
            "--suffix-output", "-o",
            type=str,
            help="Replacement suffix for output VCF files. Default: _variants_final",
            default="_variants_final"
        )
        a_parser.add_argument(
            "--dir-input", "-I",
            type=str,
            help="Directory for input VCF files. Default: current directory.",
            default="."
        )
        a_parser.add_argument(
            "--dir-output", "-O",
            type=str,
            help="Directory for output VCF files. Default: current directory.",
            default="."
        )
    else:
        a_parser.add_argument('vcf_file', help='input VCF file')
    return a_parser.parse_args()


def var_filter(
    file_vcf,
    quality_rms_min=10,
    depth_read_min=2,
    depth_read_max=10000000,
    alternatives_min=2,
    wipe_int_around=3,
    window_size=10,
    p_min_strand=1e-4,
    p_min_baseq=1e-100,
    p_min_mapq=0.0,
    p_min_distance_end=1e-4,
    p_min_hwe=1e-4,
    mxgq_min=0,
    mxsp_max=1000,
    print_filtered=False,
    output_main=sys.stdout,
    output_filtered=sys.stderr
):
    # args = parse_args()

    # Initialize VCF reader and print header
    vcf = VCF(file_vcf)
    print(vcf.raw_header, end='', file=output_main)

    # Calculate maximum distance for staging
    ol = window_size  # Window size for adjacent gaps
    ow = wipe_int_around  # Window size around gaps for SNPs
    max_dist = max(ol, ow)

    # Staging list: [score << 2 | type, flt_tag, indel_span, variant]
    staging = []

    # Process each variant
    for variant in vcf:
        # Skip non-variant sites or those with unknown reference
        if variant.REF == 'N':
            continue
        if len(variant.ALT) == 0 or variant.ALT == ['.']:
            continue

        # Determine variant type
        ref_len = len(variant.REF)
        if ref_len == 1:
            if all(len(alt) == 1 for alt in variant.ALT):
                type_ = 1  # SNP
            else:
                type_ = 3  # Indel
        else:
            if all(len(alt) == ref_len for alt in variant.ALT):
                type_ = 2  # MNP
            else:
                type_ = 3  # Indel

        # Clear staging list of out-of-range variants
        while staging and not (
                staging[0][3].CHROM == variant.CHROM
                and
                staging[0][3].POS + staging[0][2] + max_dist >= variant.POS
        ):
            var_filter_outing(
                staging_item=staging.pop(0),
                is_print=print_filtered,
                output_main=output_main,
                output_filtered=output_filtered
            )

        # Extract depth and quality metrics
        dp4 = variant.INFO.get('DP4')
        if dp4:
            dp = sum(dp4)
            dp_alt = dp4[2] + dp4[3]
        else:
            dp = variant.INFO.get('DP', -1)
            dp_alt = -1  # Assuming dp_alt is not available without DP4

        mq = variant.INFO.get('MQ', -1)

        # Apply basic filters
        flt = 0
        if dp >= 0:
            if dp < depth_read_min:
                flt = 2  # -- below minimum depth
            elif dp > depth_read_max:
                flt = 3  # -- above maximum depth
        # noinspection PyChainedComparisons
        if dp_alt >= 0 and dp_alt < alternatives_min:
            flt = 4  # -- insufficient alternate bases
        # noinspection PyChainedComparisons
        if flt == 0 and mq >= 0 and mq < quality_rms_min:
            flt = 1  # -- below minimum mapping quality

        # Apply PV4 filter (strand, baseQ, mapQ, end distance biases)
        if flt == 0:
            pv4 = variant.INFO.get('PV4')
            if pv4 and len(pv4) == 4:
                if (
                    pv4[0] < p_min_strand or pv4[1] < p_min_baseq
                    or
                    pv4[2] < p_min_mapq or pv4[3] < p_min_distance_end
                ):
                    flt = 7

        # Apply MXGQ/MXSP filter
        if flt == 0:
            mxgq = variant.INFO.get('MXGQ')
            mxsp = variant.INFO.get('MXSP')
            if (
                (mxgq is not None and mxgq < mxgq_min)
                or
                (mxsp is not None and mxsp >= mxsp_max)
            ):
                flt = 8

        # Apply HWE filter
        if flt == 0:
            g3 = variant.INFO.get('G3')
            hwe = variant.INFO.get('HWE')
            if g3 and hwe is not None and hwe < p_min_hwe:
                p = 2 * g3[0] + g3[1]
                if 0 < p < 1:
                    f = 1 - g3[1] / (p * (1 - p))
                    if f < 0:
                        flt = 9

        # Calculate score and indel span
        qual = variant.QUAL if variant.QUAL is not None else 0
        score = int(qual * 100 + (dp_alt if dp_alt >= 0 else 0))  # Cast to int
        rlen = ref_len - 1 if type_ == 3 else 0  # Indel span for indels

        # Apply indel and proximity filters
        if flt == 0:
            if type_ == 3:  # Indel
                # Filter nearby SNPs/MNPs
                for x in staging:
                    x_score, x_flt, x_indel_span, x_variant = x
                    x_type = x_score & 3
                    if (
                        x_type != 3 and x_flt == 0
                        and
                        x_variant.POS + x_indel_span + ow >= variant.POS
                    ):
                        x[1] = 5  # SNP near indel
                # Filter overlapping indels
                for x in staging:
                    x_score, x_flt, x_indel_span, x_variant = x
                    x_type = x_score & 3
                    if (
                        x_type == 3 and x_flt == 0
                        and
                        x_variant.POS + x_indel_span + ol >= variant.POS
                    ):
                        if (x_score >> 2) < score:
                            x[1] = 6  # Lower-scoring indel
                        else:
                            flt = 6  # -- current indel filtered
                            break
            elif type_ == 1:  # SNP
                # Check for adjacent indels only within the specified window
                for x in staging:
                    x_score, x_flt, x_indel_span, x_variant = x
                    x_type = x_score & 3
                    if (
                        x_type == 3 and x_flt == 0
                        and
                        x_variant.POS + x_indel_span + ow >= variant.POS
                    ):
                        flt = 5  # -- SNP near indel
                        break
            # For SNPs/MNPs, check overlapping variants
            if type_ != 3 and flt == 0:
                for x in staging:
                    x_score, x_flt, x_indel_span, x_variant = x
                    x_type = x_score & 3
                    if (
                        x_type != 3 and x_flt == 0
                        and
                        x_variant.POS + x_indel_span >= variant.POS
                    ):
                        if (x_score >> 2) < score:
                            x[1] = 8  # Lower-scoring SNP/MNP
                        else:
                            flt = 8  # -- current SNP/MNP filtered
                            break

        # Add variant to staging list
        staging.append([score << 2 | type_, flt, rlen, variant])

    # Process remaining variants in staging
    while staging:
        var_filter_outing(
            staging_item=staging.pop(0),
            is_print=print_filtered,
            output_main=output_main,
            output_filtered=output_filtered
        )


def main(
        suffix_input: str,
        suffix_output: str,
        dir_input: str | os.PathLike,
        dir_output: str | os.PathLike,
        quality_rms_min: int,
        depth_read_min: int,
        depth_read_max: int,
        alternatives_min: int,
        wipe_int_around: int,
        window_size: int,
        p_min_strand: float,
        p_min_baseq: float,
        p_min_mapq: float,
        p_min_distance_end: float,
        p_min_hwe: float,
        mxgq_min: int,
        mxsp_max: int,
        print_filtered: bool | None,
):
    print(f'Bulk var-filtering files ends with "{suffix_input}.vcf" from "{dir_input}" .')
    print(f'Output will end with "{suffix_output}.vcf" and will be written to "{dir_output}" .')
    print()
    for path_vcf_in in glob.glob(f"*{suffix_input}.vcf", root_dir=dir_input):
        print()
        print(f'Working with "{path_vcf_in}" ... ')
        path_vcf_out = Path(dir_output, f"{Path(path_vcf_in).stem}{suffix_output}.vcf")
        with path_vcf_out.open('w') as f_out:
            var_filter(
                file_vcf=path_vcf_in,
                quality_rms_min=quality_rms_min,
                depth_read_min=depth_read_min,
                depth_read_max=depth_read_max,
                alternatives_min=alternatives_min,
                wipe_int_around=wipe_int_around,
                window_size=window_size,
                p_min_strand=p_min_strand,
                p_min_baseq=p_min_baseq,
                p_min_mapq=p_min_mapq,
                p_min_distance_end=p_min_distance_end,
                p_min_hwe=p_min_hwe,
                mxgq_min=mxgq_min,
                mxsp_max=mxsp_max,
                print_filtered=print_filtered,
                output_main=f_out,
                output_filtered=sys.stdout
            )
        print(f'Wrote output to "{path_vcf_out}" .')
        print()


if __name__ == '__main__':
    args = parse_args(extended=True)
    main(
        suffix_input=args.suffix_input,
        suffix_output=args.suffix_output,
        dir_input=args.dir_input,
        dir_output=args.dir_output,
        quality_rms_min=args.quality_rms_min,
        depth_read_min=args.depth_read_min,
        depth_read_max=args.depth_read_max,
        alternatives_min=args.alternatives_min,
        wipe_int_around=args.wipe_int_around,
        window_size=args.window_size,
        p_min_strand=args.p_min_strand,
        p_min_baseq=args.p_min_baseq,
        p_min_mapq=args.p_min_mapq,
        p_min_distance_end=args.p_min_distance_end,
        p_min_hwe=args.p_min_hwe,
        mxgq_min=args.mxgq_min,
        mxsp_max=args.mxsp_max,
        print_filtered=args.print_filtered
    )
