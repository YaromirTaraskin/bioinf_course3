#!/usr/bin/env python3
"""
A Python console program to analyze VCF files.
This script processes VCF files --
analyzes variants types,
analyzes and plots DP distribution.
"""
import argparse
import os
import subprocess
from collections import defaultdict
import matplotlib.pyplot as plt
from cyvcf2 import VCF

# noinspection PyUnresolvedReferences
import seaborn as sns


def wget_download_file_if_not_exists(url, path_local: "os.PathLike[str] | str"):
    """Download a file from a URL to a local path if it doesn’t exist."""
    if not os.path.exists(path_local):
        print(f'Downloading "{url}" to "{path_local}" ...')
        subprocess.run(["wget", str(url), "-O", str(path_local)], check=True)
    else:
        print(f'Using existing file "{path_local}" .')


def tabix_extract_vcf_header_to(infile_vcf: "os.PathLike[str] | str", outfile_header: "os.PathLike[str] | str"):
    """Extract the header from a VCF file using tabix."""
    print(f'Extracting header to "{outfile_header}" ...')
    header_cmd = ["tabix", "-H", infile_vcf]
    with open(outfile_header, "w") as outfile:
        subprocess.run(header_cmd, stdout=outfile, check=True)


def add_info_definitions(header_file: "os.PathLike[str] | str", new_info_lines: "list[str]"):
    """Add new INFO field definitions to the header before the #CHROM line."""
    print("Adding missing INFO definitions...")
    with open(header_file, "r") as infile:
        header_lines = infile.readlines()

    chrom_line_index = next(i for i, line in enumerate(header_lines) if line.startswith("#CHROM"))
    header_lines = header_lines[:chrom_line_index] + new_info_lines + header_lines[chrom_line_index:]

    with open(header_file, "w") as outfile:
        outfile.writelines(header_lines)


def tabix_extract_vcf_region_data(
        infile_vcf: "os.PathLike[str] | str", region_str: str, outfile_data: "os.PathLike[str] | str"
):
    """Extract data for a specific region from the VCF file using tabix."""
    print(f'Extracting data for region "{region_str}" ...')
    data_cmd = ["tabix", infile_vcf, region_str]
    with open(outfile_data, "w") as outfile:
        subprocess.run(data_cmd, stdout=outfile, check=True)


def vcf_combine_header_and_data_to(
        infile_header: "os.PathLike[str] | str", infile_data: "os.PathLike[str] | str",
        outfile_vcf: "os.PathLike[str] | str"
):
    """Combine the header and data into a single VCF file."""
    print(f'Creating complete VCF file as "{outfile_vcf}" ...')
    with open(outfile_vcf, "w") as outfile:
        with open(infile_header, "r") as header_in:
            outfile.write(header_in.read())
        with open(infile_data, "r") as data_in:
            outfile.write(data_in.read())


def bz_compress_and_tabix_index_vcf_to(
        vcf_original: "os.PathLike[str] | str", vcf_compressed: "os.PathLike[str] | str"
):
    """Compress the VCF file using bgzip and index it with tabix."""
    print(f'Compressing output to "{vcf_compressed}" ...')
    subprocess.run(["bgzip", "-c", vcf_original], stdout=open(vcf_compressed, "wb"), check=True)
    print("Indexing VCF file...")
    subprocess.run(["tabix", "-p", "vcf", vcf_compressed], check=True)


def remove_files(files_to_remove: "list[os.PathLike[str]] | list[str]"):
    """Remove files."""
    for file in files_to_remove:
        if os.path.exists(file):
            print(f'CLEANUP: removing file "{file}" ...')
            os.remove(file)


def download_and_preprocess_vcf_data_to(
        url_base,
        vcf_filename: "os.PathLike[str] | str",
        region_string: str,
        info_to_add: "list[str]",
        vcf_out_compressed: "os.PathLike[str] | str"
):
    """Download and prepare the data"""
    # Define file paths and URLs
    url = f"{url_base}{vcf_filename}"
    url_index = f"{url}.tbi"
    vcf_index = f"{vcf_filename}.tbi"
    vcf_out = f"{vcf_out_compressed}".rsplit(sep=".", maxsplit=1)[0]
    vcf_txt_header = f"{vcf_out}_header.txt"
    vcf_txt_data = f"{vcf_out}_data.txt"

    wget_download_file_if_not_exists(url, vcf_filename)
    wget_download_file_if_not_exists(url_index, vcf_index)

    tabix_extract_vcf_header_to(vcf_filename, vcf_txt_header)
    add_info_definitions(vcf_txt_header, info_to_add)

    tabix_extract_vcf_region_data(vcf_filename, region_string, vcf_txt_data)
    vcf_combine_header_and_data_to(vcf_txt_header, vcf_txt_data, vcf_out)

    bz_compress_and_tabix_index_vcf_to(vcf_out, vcf_out_compressed)
    remove_files([vcf_txt_header, vcf_txt_data, vcf_out])


def count_variant_types(vcf_file: "os.PathLike[str] | str"):
    """
    Count the types of variants and number of alternatives for SNPs in the VCF file.

    Args:
        vcf_file (str): Path to the VCF file.

    Returns:
        tuple: (variant_types, num_alts) where variant_types is a dict of (var_type, var_subtype) counts,
               and num_alts is a dict of alternative counts for SNPs.
    """
    print("Counting the types of variants and number of alternatives for SNPs: ")
    v = VCF(vcf_file)
    variant_types = defaultdict(int)
    num_alts = defaultdict(int)

    for i, variant in enumerate(v, start=1):
        print(f"\r  processing variant {i} ... ", end='')
        variant_types[(variant.var_type, variant.var_subtype)] += 1
        if variant.var_type == 'snp':
            num_alts[len(variant.ALT)] += 1

    print("Done.")
    return variant_types, num_alts


def analyze_dp_distribution(vcf_file: "os.PathLike[str] | str"):
    """
    Analyze the depth (DP) distribution for SNP variants with one alternative.

    Args:
        vcf_file (str): Path to the VCF file.

    Returns:
        tuple: (dps, dp_dist) where
            dps is a sorted list of DP values, and
            dp_dist is a list of their frequencies.
    """
    print("Analyzing depth (DP) distribution for SNP variants with one alternative: ")
    v = VCF(vcf_file)
    sample_dp = defaultdict(int)

    for i, variant in enumerate(v, start=1):
        print(f"\r  processing variant {i} ... ", end='')
        if not variant.is_snp or len(variant.ALT) != 1:
            continue
        for dp in variant.format('DP'):
            sample_dp[dp] += 1

    print("\rSorting...", end='')
    dps = sorted(sample_dp.keys())
    dp_dist = [sample_dp[x] for x in dps]
    print("\rDP distribution analysis complete.")
    return dps, dp_dist


def plot_dp_distribution(dps, dp_dist, plot_limit=None, figname='dp_distribution.png'):
    """
    Plot the DP distribution and save it to a file.

    Args:
        dps (list): Sorted list of DP values.
        dp_dist (list): List of frequencies corresponding to dps.
        plot_limit (int): maximum values to plot from the start
        figname (str | os.PathLike): The name of the file to save the plot to.
    """
    print("Generating DP distribution plot...")
    fig, ax = plt.subplots(figsize=(16, 9))

    # Plot the first plot_limit DP values (or all if fewer)
    plot_limit = len(dp_dist) if not plot_limit else min(plot_limit, len(dp_dist))
    ax.plot(dp_dist[:plot_limit], 'r', label='DP Frequency')

    # Add a vertical line at the mode
    mode_index = dp_dist.index(max(dp_dist))
    ax.axvline(mode_index, color='b', linestyle='--', label=f'Mode DP = {dps[mode_index]}')

    # Customize the plot
    ax.set_xlabel('DP Value')
    ax.set_ylabel('Frequency')
    ax.set_title('DP Distribution for SNPs (First 50 Values)')
    ax.legend()

    # Save the plot
    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    print(f'DP distribution plot saved to "{figname}".')


def analyze_variant_types(vcf_compressed: "os.PathLike[str] | str"):
    """Analyze variant types

    Args:
        vcf_compressed (str | os.PathLike): Path to the gunzip-compressed VCF file
    """
    variant_types, num_alts = count_variant_types(vcf_compressed)
    print()
    print("Variant Types:")
    for (var_type, var_subtype), count in variant_types.items():
        print(f"  {var_type} ({var_subtype}): {count}")
    print()
    print("Number of Alternatives for SNPs:")
    for alt_count, freq in num_alts.items():
        print(f"  {alt_count} alternatives: {freq}")


def analyze_and_plot_dp_distribution(vcf_compressed: "os.PathLike[str] | str", plot_limit=None):
    """Analyze and plot DP distribution

    Args:

        vcf_compressed (str | os.PathLike): Path to the gunzip-compressed VCF file
        plot_limit (int): maximum values to plot from the start
    """
    dps, dp_dist = analyze_dp_distribution(vcf_compressed)
    mode_index = dp_dist.index(max(dp_dist))
    mode_dp = dps[mode_index]
    print()
    print(f"The DP value with the highest frequency is {mode_dp} with {max(dp_dist)} occurrences.")
    plot_dp_distribution(dps, dp_dist, plot_limit)


def main(*, url_base, vcf_filename, region, info_missed_in_original, vcf_compressed, dp_plot_limit):
    """Main function to execute the VCF analysis."""

    download_and_preprocess_vcf_data_to(
        url_base, vcf_filename, region, info_missed_in_original, vcf_compressed
    )
    if not os.path.exists(vcf_compressed):
        print(f'Error: VCF file "{vcf_compressed}" not found.')
        return

    analyze_variant_types(vcf_compressed)
    analyze_and_plot_dp_distribution(vcf_compressed, plot_limit=dp_plot_limit)


if __name__ == '__main__':
    a_parser = argparse.ArgumentParser(
        description="Analyze and plot DP distribution for SNPs in a VCF file."
    )
    a_parser.add_argument(
        '-u', '--urlfrom',
        default='ftp://ftp.1000genomes.ebi.ac.uk/'
                'vol1/ftp/release/20130502/supporting/'
                'vcf_with_sample_level_annotation/',
        help='URL from which to download the VCF file, default: '
             '"ftp://ftp.1000genomes.ebi.ac.uk/'
             'vol1/ftp/release/20130502/supporting/'
             'vcf_with_sample_level_annotation/"'
    )
    a_parser.add_argument(
        '-f', '--original',
        default='ALL.chr22.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.vcf.gz',
        help='Original VCF file name on the server, default: '
             '"ALL.chr22.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.vcf.gz"'
    )
    a_parser.add_argument(
        '-r', '--region',
        default='22:1-17000000',
        help='Region to analyze, default: "22:1-17000000"'
    )
    a_parser.add_argument(
        '-p', '--prepared',
        default='genotypes.vcf.gz',
        help='Path for the prepared VCF file to analyze, default: "genotypes.vcf.gz"'
    )
    a_parser.add_argument(
        '-l', '--limit',
        default=50,
        help='Limit for DP distribution plot, default: 50'
    )

    args = a_parser.parse_args()

    main(
        url_base=args.urlfrom,
        vcf_filename=args.original,
        region=args.region,
        info_missed_in_original=[
            '##INFO=<ID=SAS_AF,Number=1,Type=Float,Description="Allele frequency in South Asian population">\n',
            '##INFO=<ID=EAS_AF,Number=1,Type=Float,Description="Allele frequency in East Asian population">\n'
        ],
        vcf_compressed=args.prepared,
        dp_plot_limit=args.limit
    )
