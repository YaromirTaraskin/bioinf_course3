#!/usr/bin/env python3
import argparse

import working_with_fastq as worfas


def download_data(urls_file):
    """
    Download FASTQ files from a list of URLs and return their local filenames.

    Args:
        urls_file (str | os.PathLike): Path to a text file containing URLs, one per line.

    Returns:
        list: A list of downloaded filenames.
    """
    filesnames_downloaded = []
    with open(urls_file) as filestream_urls:
        for url_reads_fastq in filestream_urls:
            url_reads_fastq = url_reads_fastq.strip()
            if not url_reads_fastq:
                continue
            filename_fastq_gz = url_reads_fastq.rsplit('/', 1)[-1]
            worfas.wget_download_file_if_not_exists(url_reads_fastq, filename_fastq_gz)
            filesnames_downloaded.append(filename_fastq_gz)
    return filesnames_downloaded


def main(urls_file):
    """
    Downloads FASTQ files from a list of URLs and performs per-file quality control and analysis.

    Args:
        urls_file (str | os.PathLike): Path to a text file containing URLs, one per line.

    Returns:
        None
    """
    filenames_downloaded = download_data(urls_file)
    for filename_fastq_gz in filenames_downloaded:
        worfas.process_fastq_perform_analyses(filename_fastq_gz)


if __name__ == "__main__":
    default_urls_file = 'urls.reads.untrimmed.fastq.gz.txt'
    a_parser = argparse.ArgumentParser(
        description='Manager to orchestrate reads untrimmed FASTQ file download and processing.'
    )
    a_parser.add_argument(
        "--urls-file", "-f",
        help='Text file with URLs to use, one per line. '
             'For any file in directory having same name as in the URL, '
             f'it will be used instead of downloading. Default: "{default_urls_file}"',
        default=default_urls_file
    )

    args = a_parser.parse_args()
    main(args.urls_file)
