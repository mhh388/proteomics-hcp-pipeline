"""
src/01_parse_mztab.py
─────────────────────
Step 1 of the Proteomics HCP Pipeline.

What this script does:
  1. Parses the mzTab file from PRIDE PXD000509
  2. Extracts metadata, protein table (PRT), and PSM table
  3. Annotates proteins with UniProt gene names via API
  4. Saves clean structured data to CSV and SQLite
  5. Prints a summary of HCPs identified

Dataset: PXD000509 - E. coli HCP analysis in therapeutic protein purification
Instrument: LTQ Orbitrap Velos
Search engine: Sequest via Proteome Discoverer v1.3

Run:
  python src/01_parse_mztab.py
"""

import re
import sys
import json
import time
import sqlite3
import requests
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timezone
from loguru import logger

# ── Paths ─────────────────────────────────────────────────────
BASE_DIR     = Path(__file__).resolve().parent.parent
DATA_RAW     = BASE_DIR / "data/raw"
DATA_PROC    = BASE_DIR / "data/processed"
RESULTS_DIR  = BASE_DIR / "results"
DB_PATH      = DATA_PROC / "proteomics.db"

DATA_PROC.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Find the mzTab file automatically
MZTAB_FILES = list(DATA_RAW.glob("*.mztab"))
if not MZTAB_FILES:
    raise FileNotFoundError(
        f"No .mztab file found in {DATA_RAW}\n"
        f"Run: gunzip data/raw/*.mztab.gz"
    )
MZTAB_FILE = MZTAB_FILES[0]


# ── mzTab Parser ──────────────────────────────────────────────

def parse_mztab(filepath: Path) -> dict:
    """
    Parse an mzTab file into structured sections.

    mzTab sections:
      MTD - metadata
      PRH/PRT - protein table header/rows
      PSH/PSM - peptide-spectrum match header/rows
    """
    logger.info(f"Parsing mzTab: {filepath.name}")

    metadata = {}
    prt_rows  = []
    prt_cols  = []
    psm_rows  = []
    psm_cols  = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("COM"):
                continue

            parts = line.split("\t")
            prefix = parts[0]

            if prefix == "MTD":
                if len(parts) >= 3:
                    metadata[parts[1]] = parts[2]

            elif prefix == "PRH":
                prt_cols = parts[1:]

            elif prefix == "PRT":
                prt_rows.append(parts[1:])

            elif prefix == "PSH":
                psm_cols = parts[1:]

            elif prefix == "PSM":
                psm_rows.append(parts[1:])

    # Build DataFrames
    df_proteins = pd.DataFrame(prt_rows, columns=prt_cols) if prt_rows else pd.DataFrame()
    df_psms     = pd.DataFrame(psm_rows, columns=psm_cols) if psm_rows else pd.DataFrame()

    logger.info(f"Proteins (PRT): {len(df_proteins)} rows")
    logger.info(f"PSMs:           {len(df_psms)} rows")
    logger.info(f"Metadata keys:  {len(metadata)}")

    return {
        "metadata":  metadata,
        "proteins":  df_proteins,
        "psms":      df_psms,
    }


# ── Data Cleaning ─────────────────────────────────────────────

def clean_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and type-cast the protein table."""
    df = df.copy()

    # Rename for clarity
    rename = {
        "accession":                        "uniprot_accession",
        "best_search_engine_score[1]":      "best_score",
        "search_engine_score[1]_ms_run[1]": "score_ms_run1",
        "num_psms_ms_run[1]":               "num_psms",
        "num_peptides_distinct_ms_run[1]":  "num_peptides_distinct",
        "num_peptides_unique_ms_run[1]":    "num_peptides_unique",
        "protein_coverage":                 "protein_coverage",
        "opt_global_cv_PRIDE:0000303_Decoy_hit": "is_decoy",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    # Type conversions
    for col in ["best_score", "score_ms_run1", "protein_coverage"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ["num_psms", "num_peptides_distinct", "num_peptides_unique", "is_decoy"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    # Filter out decoys
    if "is_decoy" in df.columns:
        before = len(df)
        df = df[df["is_decoy"] == 0]
        logger.info(f"Removed {before - len(df)} decoy hits")

    return df


def clean_psms(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and type-cast the PSM table."""
    df = df.copy()

    for col in ["retention_time", "charge", "exp_mass_to_charge", "calc_mass_to_charge"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Calculate mass error (ppm) — key QC metric
    if "exp_mass_to_charge" in df.columns and "calc_mass_to_charge" in df.columns:
        df["mass_error_ppm"] = (
            (df["exp_mass_to_charge"] - df["calc_mass_to_charge"])
            / df["calc_mass_to_charge"] * 1e6
        ).round(4)
        
        df.loc[df["mass_error_ppm"].abs() > 100, "mass_error_ppm"] = np.nan
        
    # Filter extreme mass error outliers (>100 ppm = bad data, not real Orbitrap measurements)
    if "mass_error_ppm" in df.columns:
        df.loc[df["mass_error_ppm"].abs() > 100, "mass_error_ppm"] = np.nan

    # Peptide length
    if "sequence" in df.columns:
        df["peptide_length"] = df["sequence"].str.len()

    return df


# ── UniProt Annotation ────────────────────────────────────────

def annotate_with_uniprot(df: pd.DataFrame, accession_col: str = "uniprot_accession") -> pd.DataFrame:
    """
    Fetch gene names and protein names from UniProt REST API.
    This turns raw accessions like P00634 into meaningful names
    like 'phoA / Bacterial alkaline phosphatase'.
    """
    logger.info("Fetching protein annotations from UniProt...")

    accessions = df[accession_col].dropna().unique().tolist()
    annotations = {}

    for acc in accessions:
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                gene_name = ""
                if data.get("genes"):
                    gene_name = data["genes"][0].get("geneName", {}).get("value", "")
                protein_name = ""
                if data.get("proteinDescription"):
                    rec = data["proteinDescription"].get("recommendedName", {})
                    protein_name = rec.get("fullName", {}).get("value", "")
                    if not protein_name:
                        sub = data["proteinDescription"].get("submittedNames", [])
                        if sub:
                            protein_name = sub[0].get("fullName", {}).get("value", "")
                organism = data.get("organism", {}).get("scientificName", "")
                annotations[acc] = {
                    "gene_name":    gene_name,
                    "protein_name": protein_name,
                    "organism":     organism,
                }
                logger.debug(f"  {acc}: {gene_name} — {protein_name[:50]}")
            else:
                annotations[acc] = {"gene_name": "", "protein_name": "", "organism": ""}
            time.sleep(0.2)  # Be polite to UniProt API
        except Exception as e:
            logger.warning(f"UniProt lookup failed for {acc}: {e}")
            annotations[acc] = {"gene_name": "", "protein_name": "", "organism": ""}

    annot_df = pd.DataFrame.from_dict(annotations, orient="index").reset_index()
    annot_df.columns = [accession_col, "gene_name", "protein_name", "organism"]

    return df.merge(annot_df, on=accession_col, how="left")


# ── HCP Classification ────────────────────────────────────────

def classify_hcps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add HCP risk classification based on known problematic HCPs
    and PSM count. This mirrors real HCP risk assessment workflows.
    """
    # Known high-risk E. coli HCPs from literature
    HIGH_RISK_GENES = {
        "phoa",   # Bacterial alkaline phosphatase
        "ompa",   # Outer membrane protein — immunogenic
        "ompc",   # Outer membrane protein C
        "ompf",   # Outer membrane porin F
        "dnak",   # Heat shock protein
        "groel",  # GroEL chaperonin
        "lpp",    # Murein lipoprotein
        "pal",    # Peptidoglycan-associated lipoprotein
    }

    def get_risk(row):
        gene = str(row.get("gene_name", "")).lower()
        psms = row.get("num_psms", 0)
    # Gene name check FIRST — takes priority over PSM count
        if gene in HIGH_RISK_GENES:
            return "High"
        elif psms >= 5:
            return "Medium"
        elif psms >= 2:
            return "Low"
        else:
            return "Trace"

    df["hcp_risk"] = df.apply(get_risk, axis=1)
    return df


# ── Database Storage ──────────────────────────────────────────

def save_to_sqlite(proteins: pd.DataFrame, psms: pd.DataFrame, metadata: dict):
    """Save parsed data to SQLite for dashboard querying."""
    conn = sqlite3.connect(DB_PATH)

    proteins.to_sql("proteins", conn, if_exists="replace", index=False)
    psms.to_sql("psms", conn, if_exists="replace", index=False)

    # Save metadata as key-value table
    meta_df = pd.DataFrame(list(metadata.items()), columns=["key", "value"])
    meta_df.to_sql("metadata", conn, if_exists="replace", index=False)

    conn.close()
    logger.success(f"Saved to SQLite: {DB_PATH}")


# ── Summary Report ────────────────────────────────────────────

def print_summary(proteins: pd.DataFrame, psms: pd.DataFrame, metadata: dict):
    """Print a human-readable summary of the HCP analysis."""
    print(f"\n{'='*60}")
    print(f"  HCP ANALYSIS SUMMARY — PXD000509")
    print(f"{'='*60}")
    print(f"  Dataset     : {metadata.get('title', 'N/A')[:50]}")
    print(f"  Instrument  : {metadata.get('instrument[1]-name', 'N/A')}")
    print(f"  Sample      : {metadata.get('sample[1]-description', 'N/A')}")
    print(f"  Species     : E. coli")
    print(f"{'='*60}")
    print(f"  Total HCPs identified : {len(proteins)}")
    print(f"  Total PSMs            : {len(psms)}")
    print(f"  Unique peptides       : {psms['sequence'].nunique() if 'sequence' in psms.columns else 'N/A'}")
    print(f"\n── Top 10 HCPs by PSM count ──────────────────────────")

    display_cols = ["uniprot_accession", "gene_name", "protein_name",
                    "num_psms", "num_peptides_distinct", "hcp_risk"]
    display_cols = [c for c in display_cols if c in proteins.columns]

    top10 = proteins.sort_values("num_psms", ascending=False).head(10)
    print(top10[display_cols].to_string(index=False))

    if "hcp_risk" in proteins.columns:
        print(f"\n── HCP Risk Distribution ──────────────────────────────")
        risk_counts = proteins["hcp_risk"].value_counts()
        for risk, count in risk_counts.items():
            print(f"  {risk:<10} : {count} proteins")

    if "mass_error_ppm" in psms.columns:
        mae = psms["mass_error_ppm"].abs().mean()
        print(f"\n── Mass Accuracy ──────────────────────────────────────")
        print(f"  Mean mass error : {mae:.3f} ppm")
        print(f"  (< 5 ppm = excellent Orbitrap performance)")

    print(f"{'='*60}\n")


# ── Main ──────────────────────────────────────────────────────

def run_parse():
    logger.info("=" * 60)
    logger.info("  STEP 1 — PARSE mzTab + ANNOTATE HCPs")
    logger.info("=" * 60)

    # 1. Parse mzTab
    parsed = parse_mztab(MZTAB_FILE)

    # 2. Clean
    proteins = clean_proteins(parsed["proteins"])
    psms     = clean_psms(parsed["psms"])

    # 3. Annotate with UniProt
    proteins = annotate_with_uniprot(proteins)

    # 4. Classify HCP risk
    proteins = classify_hcps(proteins)

    # 5. Save to CSV
    proteins.to_csv(DATA_PROC / "proteins_annotated.csv", index=False)
    psms.to_csv(DATA_PROC / "psms_clean.csv", index=False)
    logger.success(f"Saved: data/processed/proteins_annotated.csv ({len(proteins)} proteins)")
    logger.success(f"Saved: data/processed/psms_clean.csv ({len(psms)} PSMs)")

    # 6. Save to SQLite
    save_to_sqlite(proteins, psms, parsed["metadata"])

    # 7. Save metadata
    with open(RESULTS_DIR / "metadata.json", "w") as f:
        json.dump(parsed["metadata"], f, indent=2)

    # 8. Print summary
    print_summary(proteins, psms, parsed["metadata"])

    logger.success("STEP 1 COMPLETE ✓")
    return proteins, psms, parsed["metadata"]


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO", colorize=True)
    run_parse()
