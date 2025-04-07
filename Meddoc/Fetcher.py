# src/pubmed_fetcher/fetcher.py

import os
import requests
import fitz  # PyMuPDF
from bs4 import BeautifulSoup
from Bio import Entrez
import pandas as pd
from .writer import save_to_docx

Entrez.email = "shsiung@berkeley.edu"


def fetch_metadata(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]

        title = article.get("ArticleTitle", "")
        abstract = article.get("Abstract", {}).get("AbstractText", [""])[0]
        authors = ", ".join(
            [f"{a.get('LastName', '')} {a.get('Initials', '')}" for a in article.get("AuthorList", [])]
        )
        doi = ""
        for id in article.get("ELocationID", []):
            if id.attributes.get("EIdType") == "doi":
                doi = str(id)

        return {
            "title": title,
            "abstract": abstract,
            "authors": authors,
            "doi": doi,
        }
    except Exception as e:
        print(f"Error fetching metadata for PMID {pmid}: {e}")
        return None


def fetch_full_text(doi):
    if not doi:
        return ""
    try:
        headers = {"Accept": "text/html"}
        doi_url = f"https://doi.org/{doi}"
        response = requests.get(doi_url, headers=headers, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        # Try extracting paragraphs
        paragraphs = soup.find_all("p")
        text = "\n".join(p.get_text() for p in paragraphs if p.get_text())
        return text.strip()[:5000]  # Limit length
    except Exception as e:
        print(f"Error fetching full text for DOI {doi}: {e}")
        return ""


def fetch_pdf_text(doi):
    # Basic example, needs real PDF URL handling or Unpaywall API integration
    return ""  # PDF fetching placeholder


def process_pubmed_ids(csv_file, output_dir):
    df = pd.read_csv(csv_file)
    for pmid in df["PMID"]:
        meta = fetch_metadata(pmid)
        if meta:
            html_text = fetch_full_text(meta["doi"])
            pdf_text = fetch_pdf_text(meta["doi"])
            save_to_docx(pmid, meta, html_text, pdf_text, output_dir)
