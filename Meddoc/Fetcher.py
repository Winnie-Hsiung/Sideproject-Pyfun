# src/pubmed_fetcher/fetcher.py
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
    if not doi:
        return ""

    # Try Unpaywall API
    try:
        api_url = f"https://api.unpaywall.org/v2/{doi}?email=shsiung@berkeley.edu"
        response = requests.get(api_url, timeout=10)
        response.raise_for_status()
        data = response.json()
        pdf_url = data.get("best_oa_location", {}).get("url_for_pdf")
        if pdf_url:
            pdf_response = requests.get(pdf_url, timeout=15)
            with open("temp.pdf", "wb") as f:
                f.write(pdf_response.content)
            with fitz.open("temp.pdf") as doc:
                text = "\n".join(page.get_text() for page in doc)
            os.remove("temp.pdf")
            return text.strip()[:5000]
    except Exception as e:
        print(f"Unpaywall failed for DOI {doi}: {e}")

    # Optionally fallback to Sci-Hub (not recommended/legal concerns)
    try:
        scihub_url = f"https://sci-hub.se/{doi}"
        response = requests.get(scihub_url, timeout=10)
        soup = BeautifulSoup(response.text, "html.parser")
        iframe = soup.find("iframe")
        if iframe and "src" in iframe.attrs:
            pdf_link = iframe["src"]
            if not pdf_link.startswith("http"):
                pdf_link = f"https:{pdf_link}"
            pdf_response = requests.get(pdf_link, timeout=15)
            with open("temp.pdf", "wb") as f:
                f.write(pdf_response.content)
            with fitz.open("temp.pdf") as doc:
                text = "\n".join(page.get_text() for page in doc)
            os.remove("temp.pdf")
            return text.strip()[:5000]
    except Exception as e:
        print(f"Sci-Hub failed for DOI {doi}: {e}")

    return ""


def process_pubmed_ids(csv_file, output_dir):
    df = pd.read_csv(csv_file)
    for pmid in df["PMID"]:
        meta = fetch_metadata(pmid)
        if meta:
            html_text = fetch_full_text(meta["doi"])
            pdf_text = fetch_pdf_text(meta["doi"])
            save_to_docx(pmid, meta, html_text, pdf_text, output_dir)

