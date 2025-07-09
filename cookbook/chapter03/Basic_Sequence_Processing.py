#!/usr/bin/env python3
"""
A Python console program to
fetch a CDS from a GenBank record,
save it as FASTA, and process it.
"""
from Bio import Entrez, SeqIO, SeqRecord  # pip install biopython
import argparse
import sys


def fetch_genbank_record(sequence_id: str) -> SeqRecord.SeqRecord:
    """Fetch a GenBank record from NCBI using the provided sequence ID.

    Args:
        sequence_id (str): The ID of the GenBank record to fetch.

    Returns:
        SeqRecord.SeqRecord: The fetched GenBank record.

    Raises:
        Exception: If there is an error fetching the GenBank record.
    """
    try:
        print(f"Fetching GenBank record with ID {sequence_id} ...")
        handle = Entrez.efetch(db='nucleotide', id=sequence_id, rettype='gb')
        print(f"Reading GenBank record...")
        record_fetched = SeqIO.read(handle, 'gb')
        handle.close()
        return record_fetched
    except Exception as e:
        raise Exception(f"Failed to fetch GenBank record: {e}")


def extract_cds(genbank_record: SeqRecord.SeqRecord):
    """Extract the CDS sequence from a GenBank record.

    Args:
        genbank_record (SeqRecord.SeqRecord): The GenBank record to extract the CDS from.

    Returns:
        SeqRecord.SeqRecord | None: The extracted CDS record, or None if no CDS is found.
    """
    print(f"Extracting CDS from GenBank record {genbank_record.id} ...")
    for feature in genbank_record.features:
        if feature.type == 'CDS':
            location = feature.location
            cds_seq = genbank_record.seq[location.start:location.end]
            return SeqRecord.SeqRecord(
                cds_seq, id=genbank_record.id, description='CDS only')
    return None


def write_fasta(sequence_record, outfile) -> None:
    """
    Write a sequence record to a FASTA file.

    Args:
        sequence_record (SeqRecord.SeqRecord): The sequence record to write.
        outfile (str): The path to the output FASTA file.

    Raises:
        Exception: If there is an error writing the FASTA file.
    """
    try:
        print(f"Writing sequence record to {outfile} ...")
        with open(outfile, 'w') as output_handle:
            SeqIO.write([sequence_record], output_handle, 'fasta')
    except Exception as e:
        raise Exception(f'Failed to write FASTA file: "{outfile}" : {e}')


def process_sequence(fasta_filename) -> None:
    """
    Read and process the sequence from a FASTA file.

    Args:
        fasta_filename (str): The path to the input FASTA file.

    Raises:
        Exception: If there is an error processing the FASTA file.
    """
    try:
        record = SeqIO.read(fasta_filename, 'fasta')
        sequence = record.seq

        print("--- Sequence Details ---")
        print(f"Description: {record.description}")
        print(f"First 10 bases: {sequence[:10]}")
        print(f"First 12 bases: {sequence[:12]}")
        print(f"Last 12 bases: {sequence[-12:]}")

        rna = sequence.transcribe()
        print(f"RNA: {rna}")
        protein = sequence.translate()
        print(f"Protein: {protein}")

    except Exception as e:
        raise Exception(f"Failed to process FASTA file: {e}")


def main(email, record_id, output):
    """Fetches a CDS from a GenBank record, saves it as FASTA, and process it.

    Args:
        email (str): Email address for NCBI Entrez acces
        record_id (str): Sequence ID to fetch
        output (str | os.PathLike): Output FASTA file name
    """

    # Check if email is set properly
    if email == "your@email.com":
        print("Warning: Please set a valid email address for Entrez.email and rerun.")
        print("#")
    Entrez.email = email

    # Fetch GenBank record
    try:
        gb_record = fetch_genbank_record(record_id)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Extract CDS
    cds_record = extract_cds(gb_record)
    if cds_record is None:
        print("Error: No CDS feature found in the record.", file=sys.stderr)
        sys.exit(1)

    # Write CDS to FASTA file
    try:
        write_fasta(cds_record, output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Process the FASTA file
    try:
        process_sequence(output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    # Set up command-line argument parsing
    a_parser = argparse.ArgumentParser(
        description="Fetch a CDS from a GenBank record, save it as FASTA, and process it."
    )
    a_parser.add_argument(
        '-e', '--email',
        default="your@email.com",
        help="Your email address for NCBI Entrez access"
    )
    a_parser.add_argument(
        '-i', '--id',
        default='NM_002299',
        help="Sequence ID to fetch (default: NM_002299)"
    )
    a_parser.add_argument(
        '-o', '--output',
        default='lactase_cds.fasta',
        help="Output FASTA file name (e.g., example.fasta)"
    )

    args = a_parser.parse_args()
    main(email=args.email, record_id=args.id, output=args.output)
