# cli.py

import argparse
from src.pubmed_fetcher.fetcher import process_pubmed_ids

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PubMed articles and save to DOCX.")
    parser.add_argument("csv_file", help="Path to CSV file containing a column named 'PMID'")
    parser.add_argument("output_dir", help="Directory to save the resulting Word documents")

    args = parser.parse_args()

    process_pubmed_ids(args.csv_file, args.output_dir)

# -------------

# GUI version

# cli.py
import argparse
from src.pubmed_fetcher.fetcher import process_pubmed_ids

def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed articles and save to DOCX.")
    parser.add_argument("csv_file", help="Path to CSV file containing a column named 'PMID'")
    parser.add_argument("output_dir", help="Directory to save the resulting Word documents")

    args = parser.parse_args()

    process_pubmed_ids(args.csv_file, args.output_dir)

if __name__ == "__main__":
    main()
