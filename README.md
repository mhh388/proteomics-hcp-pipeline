# Proteomics HCP Analytics Pipeline

End-to-end LC-MS/MS host cell protein (HCP) analysis pipeline built on real biopharma data.

## Dataset
PXD000509 — E. coli HCP identification in therapeutic protein purification
Bomans et al., Roche Diagnostics GmbH | PRIDE Archive, EBI

## Pipeline
- Step 1: src/01_parse_mztab.py — Parse mzTab, UniProt annotation, HCP risk classification
- Step 2: src/02_hcp_analysis.py — Spectral counting quantification, mass accuracy QC, S3 upload
- Step 3: src/03_dashboard.py — Interactive Streamlit dashboard

## Key Findings
- 54 E. coli HCPs identified across 220 PSMs
- phoA (Bacterial alkaline phosphatase) confirmed as key high-risk HCP
- 99.49% PSMs within 5 ppm — excellent Orbitrap mass accuracy
- Mean mass error: -1.24 ppm

## Tech Stack
Python, pyteomics, pandas, Streamlit, Plotly, SQLite, AWS S3, UniProt REST API

## Run Locally
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python src/01_parse_mztab.py
python src/02_hcp_analysis.py --no-s3
streamlit run src/03_dashboard.py