#!/usr/bin/env python3

"""
A Python console program to analyze VCF files and generate plots.
This script processes
bi-allelic SNPs,
MQ0 distributions,
heterozygosity vs. read depth,
and variant effect vs. depth,
saving the results as PNG files.

# Step 1: Download the full VCF file
wget ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/preview/ag1000g.AC.phase1.AR1.vcf.gz
-O ag1000g.AC.phase1.AR1.vcf.gz

# Step 2: Index the VCF file or download the index
tabix -p vcf ag1000g.AC.phase1.AR1.vcf.gz
  or
wget ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/preview/ag1000g.AC.phase1.AR1.vcf.gz.tbi

# Step 3: Extract regions
tabix -h ag1000g.AC.phase1.AR1.vcf.gz 3L:1-200000 | bgzip -c > centro.vcf.gz
tabix -h ag1000g.AC.phase1.AR1.vcf.gz 3L:21000000-21200000 | bgzip -c > standard.vcf.gz

# Step 4: Index the extracted files (optional as it added to the program)
tabix -p vcf centro.vcf.gz
tabix -p vcf standard.vcf.gz

# Step 5: Run this program
"""

# --- Imports ---

import argparse
import subprocess
from collections import defaultdict
import functools
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from cyvcf2 import VCF


# --- Helper Functions ---

def tabix_index(vcf_file):
    """Index a VCF file using tabix.

    Args:
        vcf_file (str | os.PathLike): Path to the VCF file
    """
    print(f'Indexing VCF file {vcf_file} ...')
    try:
        subprocess.run(["tabix", "-f", "-p", "vcf", vcf_file], check=True)
    except Exception as exn:
        print(f'  An exception occurred while indexing: {exn}')


def group_records_into_windows(vcf_recs, win_size, fun):
    """
    Group VCF records into windows of specified size and apply a function to each record.

    Args:
        vcf_recs: VCF reader object
        win_size (int): Window size in base pairs
        fun (callable): Function to apply to each record

    Returns:
        list: List of lists, where each inner list contains results for a window
    """
    print(f"Grouping records into windows (of size {win_size}):")
    start = None
    win_res = []
    for i, rec in enumerate(vcf_recs, start=1):
        print(f"\r  processing record {i} ... ", end='')
        if not rec.is_snp or len(rec.ALT) > 1:
            continue
        if start is None:
            start = rec.POS
        my_win = 1 + (rec.POS - start) // win_size
        while len(win_res) < my_win:
            win_res.append([])
        win_res[my_win - 1].extend(fun(rec))
    print("Done.")
    return win_res


def apply_functions_to_windows(windows, functions):
    """
    Apply multiple functions to each window's data.

    Args:
        windows (list): List of window data
        functions (dict): Dictionary of function names to functions

    Returns:
        list: List of dictionaries with function results per window
    """
    fun_results = []
    for win in windows:
        my_funs = {}
        for f_name, fun in functions.items():
            if len(win) > 0:  # Check if the window has data
                # noinspection PyBroadException
                try:
                    my_funs[f_name] = fun(win)
                except Exception:
                    my_funs[f_name] = None  # Handle any computation errors
            else:
                my_funs[f_name] = None  # Assign a default value for empty windows
        fun_results.append(my_funs)
    return fun_results


def extract_sample_data(vcf_record, annot, np_int_type):
    """
    Extract sample-specific data from a record, filtering out missing values.

    Args:
        vcf_record: VCF record
        annot (str): Annotation name (e.g., 'MQ0')
        np_int_type: NumPy integer type for missing value comparison

    Returns:
        list: List of values for the annotation across samples
    """
    return [v for v in vcf_record.format(annot) if v > np.iinfo(np_int_type).min]


def compute_sample_attribute_relation(vcf_recs, f_1, f_2):
    """
    Compute the relationship between two sample-specific attributes.

    Args:
        vcf_recs: VCF reader object
        f_1 (callable): Function to compute first attribute (e.g., heterozygosity)
        f_2 (callable): Function to compute second attribute (e.g., DP)

    Returns:
        defaultdict: Counts of (f1, f2) pairs
    """
    print("Computing sample attribute relation:")
    rel = defaultdict(int)
    for i, rec in enumerate(vcf_recs, start=1):
        print(f"\r  processing record {i} ... ", end='')
        if not rec.is_snp:
            continue
        for pos in range(len(rec.genotypes)):
            val_1 = f_1(rec, pos)
            val_2 = f_2(rec, pos)
            if val_1 is None or val_2 == np.iinfo(type(val_2)).min:
                continue
            rel[(val_1, val_2)] += 1
    print("Done.")
    return rel


def compute_variant_attribute_relation(vcf_recs, f_1, f_2):
    """
    Compute the relationship between two variant-specific attributes.

    Args:
        vcf_recs: VCF reader object
        f_1 (callable): Function to compute first attribute (e.g., effect type)
        f_2 (callable): Function to compute second attribute (e.g., DP)

    Returns:
        defaultdict: Counts of (f1, f2) pairs
    """
    print("Computing variant attribute relation:")
    rel = defaultdict(int)
    for i, rec in enumerate(vcf_recs, start=1):
        print(f"\r  processing record {i} ... ", end='')
        if not rec.is_snp:
            continue
        # noinspection PyBroadException
        try:
            val_1 = f_1(rec)
            val_2 = f_2(rec)
            if val_1 is None or val_2 is None:
                continue
            rel[(val_1, val_2)] += 1
        except Exception:
            pass
    print("Done.")
    return rel


def convert_effect_to_index(rec, accepted_effects):
    """
    Convert variant effect annotation to an integer index.

    Args:
        rec: VCF record
        accepted_effects (list): List of accepted effect types

    Returns:
        int: Index in accepted_effects or len(accepted_effects) for 'OTHER'
    """
    try:
        annot = rec.INFO['EFF']
        master_type = annot.split('(')[0]
        return accepted_effects.index(master_type)
    except ValueError:
        return len(accepted_effects)


# --- Main Step Functions ---

def plot_snp_counts(windows, window_size, figname='SNPs.png'):
    """Plot the number of bi-allelic SNPs per window for each VCF file.

    Args:
        windows (dict[str, list]): Dictionary of VCF file names to window data
        window_size (int): Size of the window in base pairs
        figname (str | os.PathLike): Name of the output plot to save as image file
    """
    print('Plotting bi-allelic SNPs per window...', end='')
    stats = {}
    fig, ax = plt.subplots(figsize=(16, 9), dpi=300, tight_layout=True)
    for w_name, wins_data in windows.items():
        stats[w_name] = apply_functions_to_windows(wins_data, {'sum': sum})
        x_lim = [i * window_size for i in range(len(stats[w_name]))]
        ax.plot(x_lim, [x['sum'] for x in stats[w_name]], label=w_name)
    ax.legend()
    ax.set_xlabel('Genomic location in the downloaded segment', fontsize='xx-large')
    ax.set_ylabel('Number of variant sites (bi-allelic SNPs)', fontsize='xx-large')
    fig.suptitle('Number of bi-allelic SNPs along the genome', fontsize='xx-large')
    fig.savefig(figname)
    print(f'\rSaved plot of bi-allelic SNPs to "{figname}" .')
    plt.close(fig)


def process_and_plot_snp_counts(vcf_files, window_size, figname='SNPs.png'):
    """
    Plot the number of bi-allelic SNPs per window for each VCF file.

    Args:
        vcf_files (list): List of VCF file names
        window_size (int): Size of the window in base pairs
        figname (str | os.PathLike): Name of the output plot to save as image file
    """
    print()
    print(f"Processing SNP counts per window of size {window_size}")
    print()
    windows = {
        vcf_f: group_records_into_windows(
            VCF(vcf_f), window_size, lambda x: [1]
        ) for vcf_f in vcf_files
    }

    plot_snp_counts(windows, window_size, figname)
    print("--------")
    return windows


def plot_mq0_distribution(mq0_wins, window_size, figname='MQ0.png'):
    """Plot the MQ0 distribution per window

    Args:
        mq0_wins (dict[str, list]): Dictionary of VCF file names to window data
        window_size (int): Size of the window in base pairs
        figname (str | os.PathLike): Name of the output plot to save as image file
    """
    print("Plotting MQ0 distribution per window... ", end='')
    stats = {}
    colors = ['b', 'g']
    fig, ax = plt.subplots(figsize=(16, 9), dpi=300, tight_layout=True)
    for i, (w_name, wins_data) in enumerate(mq0_wins.items()):
        stats[w_name] = apply_functions_to_windows(wins_data, {
            'median': np.median,
            '95': functools.partial(np.percentile, q=95)
        })
        x_lim = [j * window_size for j in range(len(stats[w_name]))]
        ax.plot(x_lim, [x['median'] for x in stats[w_name]], label=w_name, color=colors[i])
        ax.plot(x_lim, [x['95'] for x in stats[w_name]], '--', color=colors[i])
    ax.legend()
    ax.set_xlabel('Genomic location in the downloaded segment', fontsize='xx-large')
    ax.set_ylabel('MQ0', fontsize='xx-large')
    fig.suptitle('Distribution of MQ0 along the genome', fontsize='xx-large')
    fig.savefig(figname)
    print(f'\rSaved plot of MQ0 distribution to "{figname}" .')
    plt.close(fig)


def process_and_plot_mq0_distribution(vcf_files, window_size, figname='MQ0.png'):
    """
    Plot the distribution of MQ0 per window for each VCF file.

    Args:
        vcf_files (list): List of VCF file names
        window_size (int): Size of the window in base pairs
        figname (str | os.PathLike): Name of the output plot to save as image file
    """
    print()
    print(f"Processing MQ0 distribution per window of size {window_size}")
    print()
    mq0_wins = {
        vcf_f: group_records_into_windows(
            VCF(vcf_f), window_size,
            functools.partial(
                extract_sample_data, annot='MQ0', np_int_type=np.int32
            )
        ) for vcf_f in vcf_files
    }

    plot_mq0_distribution(mq0_wins, window_size, figname)
    print("--------")
    return mq0_wins


def plot_heterozygosity_vs_depth(rels, figname='hz_vs_dp.png'):
    """
    Plot the number of calls per depth and the fraction of calls which are heterozygote (Hz).

    Args:
        rels (dict): A dictionary where the keys are file names and the values are dictionaries.
            The inner dictionaries have keys which are tuples of (heterozygote, depth) and values
            which are counts of calls with that heterozygote and depth.
        figname (str | os.PathLike): The name of the file to save the plot to.
    """
    print("Plotting heterozygosity vs depth...", end='')
    fig, ax = plt.subplots(figsize=(16, 9), dpi=300, tight_layout=True)
    ax2 = ax.twinx()
    for filename, relation in rels.items():
        dps = sorted(set(x[1] for x in relation.keys()))
        frac_hetzyg = []
        count_dp = []
        for dp in dps:
            hetzyg = 0.0
            count = 0
            for hetzyg_key, dp_key in relation.keys():
                if dp_key != dp:
                    continue
                count += relation[(hetzyg_key, dp)]
                if hetzyg_key == 1:
                    hetzyg += relation[(hetzyg_key, dp)]
            frac_hetzyg.append(hetzyg / count if count > 0 else 0)
            count_dp.append(count)
        ax.plot(dps, frac_hetzyg, label=filename)
        ax2.plot(dps, count_dp, '--', label=filename)
    ax.set_xlim(0, 75)
    ax.set_ylim(0, 0.2)
    ax2.set_ylabel('Quantity of calls', fontsize='xx-large')
    ax.set_ylabel('Fraction of Heterozygote calls', fontsize='xx-large')
    ax.set_xlabel('Sample Read Depth (DP)', fontsize='xx-large')
    ax.legend()
    fig.suptitle('Number of calls per depth and fraction of calls which are Hz', fontsize='xx-large')
    fig.savefig(figname)
    print(f'\rSaved plot of heterozygosity vs depth to "{figname}" .')
    plt.close(fig)


def process_and_plot_heterozygosity_vs_depth(vcf_files, figname='hz_vs_dp.png'):
    """
    Plot the fraction of heterozygote calls and number of calls per sample DP for each VCF file.

    Args:
        vcf_files (list): List of VCF file names
        figname (str | os.PathLike): Name of the output figure
    """
    print()
    print("Processing sample relation (heterozygosity vs. DP)")
    print()
    rels = {
        vcf_f: compute_sample_attribute_relation(
            VCF(vcf_f),
            lambda rec, pos: 1 if rec.genotypes[pos][0] != rec.genotypes[pos][1] else 0,
            lambda rec, pos: rec.format('DP')[pos][0]
        ) for vcf_f in vcf_files
    }

    plot_heterozygosity_vs_depth(rels, figname)
    print("--------")
    return rels


def plot_effect_vs_depth(eff_dp_rel, accepted_effects, figname='eff_vs_dp.png'):
    """Plot the distribution of variant DP per SNP type.

    Args:
        eff_dp_rel (dict): A dictionary where keys are tuples (effect, depth) and values are counts.
        accepted_effects (list): A list of accepted effect types.
        figname (str | os.PathLike): The name of the output figure file.
    """
    print("Plotting effect vs depth... ", end='')
    fig, ax = plt.subplots(figsize=(16, 9), dpi=300, tight_layout=True)
    boxplot_vals = [[] for _ in range(len(accepted_effects) + 1)]
    for eff_and_dp, count in eff_dp_rel.items():  # Assuming eff_dp_rel is your data
        eff, dp = eff_and_dp
        boxplot_vals[eff].extend([dp] * count)
    sns.boxplot(data=boxplot_vals, showfliers=False, ax=ax)

    # Fix the ticks and labels
    num_categories = len(accepted_effects) + 1
    ax.set_xticks(range(num_categories))  # Set tick positions
    ax.set_xticklabels(accepted_effects + ['OTHER'])  # Set labels
    ax.set_ylabel('DP (variant)', fontsize='xx-large')
    fig.suptitle('Distribution of variant DP per SNP type', fontsize='xx-large')
    fig.savefig(figname)
    print(f'\rSaved plot of effect vs depth to "{figname}" .')
    plt.close(fig)  # Close the figure to free memory


def process_and_plot_effect_vs_depth(vcf_file, accepted_effects, figname='eff_vs_dp.png'):
    """
    Plot the distribution of variant DP per SNP type for the given VCF file.

    Args:
        vcf_file (str): Name of the VCF file
        accepted_effects (list): List of accepted effect types
        figname: Name of the output figure
    """
    print()
    print(f'Processing variant relation (effect vs. DP) for {vcf_file}')
    print()
    eff_dp_rel = compute_variant_attribute_relation(
        VCF(vcf_file),
        lambda rec: convert_effect_to_index(rec, accepted_effects),
        lambda rec: int(rec.INFO['DP'])
    )

    plot_effect_vs_depth(eff_dp_rel, accepted_effects, figname)

    print("--------")
    return eff_dp_rel


# --- Main Function ---

def main(
    file_centro, file_standard, window_size_snps, window_size_mq0
):
    """
    Main function to orchestrate VCF file analysis and plot generation.

    Args:
        file_centro (str): Name of the VCF file for the 'centro' region
        file_standard (str): Name of the VCF file for the 'standard' region
        window_size_snps (int): Window size for SNP count analysis
        window_size_mq0 (int): Window size for MQ0 distribution analysis
    """
    vcf_files = [file_centro, file_standard]
    accepted_effects = ['INTERGENIC', 'INTRON', 'NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING']

    matplotlib.use('Agg')  # Use non-interactive backend

    tabix_index(file_centro)
    tabix_index(file_standard)
    print()
    print("========================================")
    print()
    print("Starting VCF file analysis...")
    print("--------")

    process_and_plot_snp_counts(vcf_files, window_size_snps, figname='bi.png')
    process_and_plot_mq0_distribution(vcf_files, window_size_mq0, figname='MQ0.png')

    process_and_plot_heterozygosity_vs_depth(vcf_files, figname='hz.png')
    process_and_plot_effect_vs_depth(file_standard, accepted_effects, figname='eff.png')

    print()
    print("Analysis complete. All plots have been saved.")


if __name__ == "__main__":
    a_parser = argparse.ArgumentParser(
        description="Processes "
                    "bi-allelic SNPs, "
                    "MQ0 distributions, "
                    "heterozygosity vs. read depth, "
                    "and variant effect vs. depth (only for 'standard' region), "
                    "saving the results as PNG files."
    )
    a_parser.add_argument(
        '-c', '--file-centro',
        default='centro.vcf.gz',
        help="Path of the VCF file for the 'centro' region, default: centro.vcf.gz"
    )
    a_parser.add_argument(
        '-f', '--file-standard',
        default='standard.vcf.gz',
        help="Path of the VCF file for the 'standard' region, default: standard.vcf.gz"
    )
    a_parser.add_argument(
        '-s', '--winsize-snps',
        default=2000,
        help="Window size for bi-allelic SNPs analysis, default: 2000"
    )
    a_parser.add_argument(
        '-m', '--winsize-mq0',
        default=5000,
        help="Window size for MQ0 distribution analysis, default: 5000"
    )
    args = a_parser.parse_args()

    main(
        file_centro=args.file_centro,
        file_standard=args.file_standard,
        window_size_snps=args.winsize_snps,
        window_size_mq0=args.winsize_mq0
    )
