#!/usr/bin/env python3
"""
A Python console program to fetch record from the NCBI databases and display information about them.
"""
import argparse
from Bio import Entrez, Medline, SeqIO, SeqRecord  # pip install biopython
import sys


def search_nucleotide(search_term: str, database: str = "nucleotide") -> list[str]:
    """Search the specified database with the given term and return a list of IDs.

    Args:
        search_term (str): The search term to query the database.
        database (str): The database to search in. Defaults to "nucleotide".

    Returns:
        list[str]: A list of record IDs found by the search.
    """
    try:
        print(f"Searching {database} for '{search_term}'...")
        handle = Entrez.esearch(db=database, term=search_term)
        rec_list = Entrez.read(handle)
        count = int(rec_list['Count'])
        print(f"Found {count} records.")
        if count > int(rec_list['RetMax']):
            print("More records available, adjusting retmax...")
            handle = Entrez.esearch(db=database, term=search_term, retmax=count)
            rec_list = Entrez.read(handle)
        return rec_list['IdList']
    except Exception as e:
        print(f"Error during search: {e}")
        return []


def fetch_records(
        id_list: list[str], database: str = "nucleotide", rettype: str = "gb"
) -> "list[SeqRecord.SeqRecord]":
    """Fetch records from the specified database using the provided ID list.

    Args:
        id_list: List of record IDs to fetch.
        database: The database to fetch records from. Defaults to "nucleotide".
        rettype: The return type of the records. Defaults to "gb".

    Returns:
        list[SeqRecord.SeqRecord]: The fetched records.
    """
    try:
        print(f"Fetching {len(id_list)} records from {database}... ")
        handle = Entrez.efetch(
            db=database, id=id_list, rettype=rettype, retmax=len(id_list)
        )
        print("  parsing (it's often takes much time) ... ")
        records_fetched = list(SeqIO.parse(handle, 'gb'))
        print(f"{len(records_fetched)} records fetched successfully.")
        return records_fetched
    except Exception as e:
        print(f"Error fetching records: {e}")
        return []


def fetch_pubmed_details(id_pubmed: str) -> None:
    """Fetch and display PubMed details for a given PubMed ID.

    Args:
        id_pubmed (str): The PubMed ID to fetch details for.
    """
    try:
        print("Fetching PubMed details... ", end='')
        handle_pubmed_medline = (
            Entrez.efetch(
                db="pubmed", id=id_pubmed, rettype="medline", retmode="text"
            )
        )
        print(' parsing... ', end='')
        medline_records = Medline.parse(handle_pubmed_medline)
        print("\rPubMed details:")
        for med_rec in medline_records:
            for rec_name, rec_value in med_rec.items():
                print(f"  {rec_name}: {rec_value}")
    except Exception as e:
        print(f"Error fetching PubMed details: {e}")


def process_record(record_to_process):
    """
    Process and display details of a single SeqRecord.

    Args:
        record_to_process (SeqRecord.SeqRecord): The SeqRecord to process.
    """
    print()
    print(f"Processing record: {record_to_process.name}")
    print(f"Description: {record_to_process.description}")

    print()
    print("Features:")
    for record_feature in record_to_process.features:
        if record_feature.type == 'gene':
            gene_name = record_feature.qualifiers.get('gene', ['Unknown'])[0]
            print(f"  Gene: {gene_name}")
        elif record_feature.type == 'exon':
            exon_loc = record_feature.location
            print(
                f"  Exon: {exon_loc.start} - {exon_loc.end}, "
                f"strand: {exon_loc.strand}"
            )
        else:
            print(f"  Other feature (not processed): {record_feature.type}")

    print()
    print("Annotations:")
    for annotation_name, annotation_value in record_to_process.annotations.items():
        print(f"  {annotation_name}: {annotation_value}")

    print()
    print(f"Sequence length: {len(record_to_process.seq)}")

    print()
    print("References:")
    print("-----------")
    references = record_to_process.annotations.get('references', [])
    for ref in references:
        print(f"Title: {ref.title}")
        print(f"Authors: {ref.authors}")
        print(f"Journal: {ref.journal}")
        if ref.pubmed_id:
            print(f"PubMed ID: {ref.pubmed_id}")
            fetch_pubmed_details(ref.pubmed_id)
        print("---")


def main(email, search_term, target_name):
    """Main function to orchestrate the fetching and processing of genetic data."""

    Entrez.email = email
    # Check if email is set properly
    if Entrez.email == "your_email@example.com":
        print("Warning: Please set a valid email address for Entrez.email and rerun.")
        print("#")

    # Step 1: Search for records
    id_list = search_nucleotide(search_term)
    if not id_list:
        print("No records found. Exiting.")
        sys.exit(1)

    # Step 2: Fetch the records
    records = fetch_records(id_list)
    if not records:
        print("Failed to fetch records. Exiting.")
        sys.exit(1)

    # Step 3: Find and process the specific record
    record_found = False
    for rec in records:
        if rec.name == target_name:
            process_record(rec)
            record_found = True
    if not record_found:
        print(f"Record with name {target_name} not found.")


if __name__ == "__main__":
    a_parser = argparse.ArgumentParser(
        description="Fetch and process genetic data from NCBI databases."
    )
    a_parser.add_argument(
        '-e', '--email',
        default="your_email@example.com",
        help="Your email address for NCBI Entrez access (Entrez.email), it's better to specify it."
    )
    a_parser.add_argument(
        '-s', '--searchterm',
        default='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]',
        help='Search term for NCBI Entrez, default: CRT[Gene Name] AND \"Plasmodium falciparum\"[Organism]'
    )
    a_parser.add_argument(
        '-t', '--target',
        default='KM288867',
        help='Target name for NCBI Entrez, default: KM288867'
    )

    args = a_parser.parse_args()

    main(
        email=args.email,
        search_term=args.searchterm,
        target_name=args.target
    )
