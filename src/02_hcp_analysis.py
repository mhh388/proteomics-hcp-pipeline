"""
src/02_hcp_analysis.py
──────────────────────
Step 2 of the Proteomics HCP Pipeline.

What this script does:
  1. Loads parsed protein and PSM data
  2. Calculates HCP abundance metrics (spectral counting)
  3. Performs statistical analysis on mass accuracy
  4. Identifies co-eluting peptides and charge state distributions
  5. Generates a structured analysis report
  6. Uploads results to S3 (optional)

Spectral counting is used as a semi-quantitative method here —
more PSMs = higher protein abundance. This mirrors real HCP
quantification workflows in biopharma.

Run:
  python src/02_hcp_analysis.py
  python src/02_hcp_analysis.py --no-s3
"""

import sys
import json
import boto3
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timezone
from loguru import logger
from dotenv import load_dotenv

BASE_DIR = Path(__file__).resolve().parent.parent
load_dotenv(BASE_DIR / ".env")

DATA_PROC   = BASE_DIR / "data/processed"
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

import os
S3_BUCKET             = os.getenv("S3_BUCKET", "healthcare-etl-mengqi")
S3_PROTEOMICS_PREFIX  = os.getenv("S3_PROTEOMICS_PREFIX", "proteomics/PXD000509")
AWS_REGION            = os.getenv("AWS_REGION", "us-east-1")


# ── Spectral Counting Quantification ─────────────────────────

def calculate_spectral_counts(proteins: pd.DataFrame, psms: pd.DataFrame) -> pd.DataFrame:
    """
    Spectral counting: use PSM count as proxy for protein abundance.
    This is the standard semi-quantitative approach for HCP analysis
    when no stable isotope labeling is used.
    """
    logger.info("Calculating spectral counting metrics...")

    # PSM-based abundance
    psm_counts = psms.groupby("accession").agg(
        total_psms         = ("sequence", "count"),
        unique_peptides    = ("sequence", "nunique"),
        avg_charge         = ("charge", "mean"),
        avg_mass_error_ppm = ("mass_error_ppm", lambda x: x.abs().mean()
                              if "mass_error_ppm" in psms.columns else np.nan),
        min_peptide_len    = ("peptide_length", "min"),
        max_peptide_len    = ("peptide_length", "max"),
        avg_peptide_len    = ("peptide_length", "mean"),
    ).reset_index()

    # Merge with protein table
    df = proteins.merge(
        psm_counts,
        left_on="uniprot_accession",
        right_on="accession",
        how="left"
    ).drop(columns=["accession"], errors="ignore")

    # Normalized spectral abundance factor (NSAF)
    # NSAF = (PSMs/protein_length) / sum(PSMs/protein_length)
    # We don't have protein length from mzTab so use raw PSM fraction
    total_psms = df["total_psms"].sum()
    df["spectral_abundance_pct"] = (df["total_psms"] / total_psms * 100).round(3)

    # Rank by abundance
    df["abundance_rank"] = df["total_psms"].rank(ascending=False, method="min").astype(int)

    logger.info(f"Top HCP: {df.sort_values('total_psms', ascending=False).iloc[0]['gene_name']} "
                f"({df.sort_values('total_psms', ascending=False).iloc[0]['total_psms']} PSMs)")

    return df.sort_values("total_psms", ascending=False).reset_index(drop=True)


def analyze_mass_accuracy(psms: pd.DataFrame) -> dict:
    """
    Analyze mass measurement accuracy — a key data quality metric
    for Orbitrap instruments. < 5 ppm is excellent.
    """
    if "mass_error_ppm" not in psms.columns:
        return {}

    errors = psms["mass_error_ppm"].dropna()
    return {
        "mean_ppm":    round(float(errors.mean()), 4),
        "median_ppm":  round(float(errors.median()), 4),
        "std_ppm":     round(float(errors.std()), 4),
        "within_5ppm": round(float((errors.abs() <= 5).mean() * 100), 2),
        "within_10ppm": round(float((errors.abs() <= 10).mean() * 100), 2),
        "min_ppm":     round(float(errors.min()), 4),
        "max_ppm":     round(float(errors.max()), 4),
    }


def analyze_charge_states(psms: pd.DataFrame) -> dict:
    """Charge state distribution — reflects peptide size and digestion efficiency."""
    if "charge" not in psms.columns:
        return {}
    charge_dist = psms["charge"].value_counts().sort_index()
    return {str(int(k)): int(v) for k, v in charge_dist.items() if pd.notna(k)}


def analyze_peptide_lengths(psms: pd.DataFrame) -> dict:
    """Peptide length distribution — quality check for tryptic digestion."""
    if "peptide_length" not in psms.columns:
        return {}
    lengths = psms["peptide_length"].dropna()
    return {
        "mean_length":    round(float(lengths.mean()), 2),
        "median_length":  round(float(lengths.median()), 2),
        "min_length":     int(lengths.min()),
        "max_length":     int(lengths.max()),
        "pct_7_to_25":   round(float(((lengths >= 7) & (lengths <= 25)).mean() * 100), 2),
    }


# ── HCP-Specific Analysis ─────────────────────────────────────

def identify_key_hcps(df: pd.DataFrame) -> dict:
    """
    Flag proteins of specific biological/regulatory significance.
    Based on published HCP risk assessment frameworks.
    """
    findings = {}

    # Bacterial alkaline phosphatase — key finding of this paper
    bap = df[df["gene_name"].str.lower() == "phoa"] if "gene_name" in df.columns else pd.DataFrame()
    if not bap.empty:
        findings["bacterial_alkaline_phosphatase"] = {
            "accession":     bap.iloc[0].get("uniprot_accession", "P00634"),
            "gene":          "phoA",
            "psms":          int(bap.iloc[0].get("total_psms", 0)),
            "rank":          int(bap.iloc[0].get("abundance_rank", 0)),
            "significance":  "Most abundant HCP — key finding. Enzymatic activity monitored by cobas assay.",
            "risk":          "High",
        }

    # Top 5 most abundant HCPs
    top5 = df.head(5)[["uniprot_accession", "gene_name", "protein_name",
                        "total_psms", "spectral_abundance_pct", "hcp_risk"]].to_dict("records")
    findings["top_5_abundant_hcps"] = top5

    # Risk distribution
    if "hcp_risk" in df.columns:
        findings["risk_distribution"] = df["hcp_risk"].value_counts().to_dict()

    return findings


# ── S3 Upload ─────────────────────────────────────────────────

def upload_results_to_s3(results_dir: Path):
    """Upload processed results to S3 proteomics data lake."""
    try:
        s3 = boto3.client("s3", region_name=AWS_REGION)
        uploaded = []
        for f in results_dir.glob("*.json"):
            key = f"{S3_PROTEOMICS_PREFIX}/results/{f.name}"
            s3.upload_file(str(f), S3_BUCKET, key)
            uploaded.append(key)
            logger.info(f"  → s3://{S3_BUCKET}/{key}")
        for f in (BASE_DIR / "data/processed").glob("*.csv"):
            key = f"{S3_PROTEOMICS_PREFIX}/processed/{f.name}"
            s3.upload_file(str(f), S3_BUCKET, key)
            uploaded.append(key)
        logger.success(f"Uploaded {len(uploaded)} files to S3")
    except Exception as e:
        logger.warning(f"S3 upload skipped: {e}")


# ── Main ──────────────────────────────────────────────────────

def run_analysis(upload_s3: bool = True):
    logger.info("=" * 60)
    logger.info("  STEP 2 — HCP QUANTIFICATION & ANALYSIS")
    logger.info("=" * 60)

    # Load parsed data
    proteins_path = DATA_PROC / "proteins_annotated.csv"
    psms_path     = DATA_PROC / "psms_clean.csv"

    if not proteins_path.exists():
        raise FileNotFoundError("Run 01_parse_mztab.py first")

    proteins = pd.read_csv(proteins_path)
    psms     = pd.read_csv(psms_path)
    logger.info(f"Loaded {len(proteins)} proteins, {len(psms)} PSMs")

    # Analyses
    df_quant      = calculate_spectral_counts(proteins, psms)
    mass_accuracy = analyze_mass_accuracy(psms)
    charge_dist   = analyze_charge_states(psms)
    peptide_stats = analyze_peptide_lengths(psms)
    key_hcps      = identify_key_hcps(df_quant)

    # Save quantified proteins
    df_quant.to_csv(DATA_PROC / "proteins_quantified.csv", index=False)

    # Save full analysis report
    report = {
        "project":         "PXD000509",
        "dataset_title":   "HCP identification in therapeutic protein purification (E. coli)",
        "instrument":      "LTQ Orbitrap Velos",
        "search_engine":   "Sequest / Proteome Discoverer v1.3",
        "total_proteins":  len(df_quant),
        "total_psms":      len(psms),
        "unique_peptides": int(psms["sequence"].nunique()) if "sequence" in psms.columns else 0,
        "mass_accuracy":   mass_accuracy,
        "charge_distribution": charge_dist,
        "peptide_length_stats": peptide_stats,
        "key_hcp_findings": key_hcps,
        "analyzed_at":     datetime.now(timezone.utc).isoformat(),
    }

    report_path = RESULTS_DIR / "hcp_analysis_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    # Print key findings
    print(f"\n{'='*60}")
    print(f"  HCP ANALYSIS RESULTS — PXD000509")
    print(f"{'='*60}")
    print(f"  Proteins identified : {report['total_proteins']}")
    print(f"  Total PSMs          : {report['total_psms']}")
    print(f"  Unique peptides     : {report['unique_peptides']}")
    if mass_accuracy:
        print(f"  Mean mass error     : {mass_accuracy['mean_ppm']} ppm")
        print(f"  PSMs within 5 ppm   : {mass_accuracy['within_5ppm']}%")
    print(f"\n── Top 5 Most Abundant HCPs ───────────────────────────")
    cols = ["uniprot_accession", "gene_name", "total_psms",
            "spectral_abundance_pct", "hcp_risk"]
    cols = [c for c in cols if c in df_quant.columns]
    print(df_quant[cols].head(5).to_string(index=False))
    print(f"\n  Full report: {report_path}")
    print(f"{'='*60}\n")

    # Upload to S3
    if upload_s3:
        upload_results_to_s3(RESULTS_DIR)

    logger.success("STEP 2 COMPLETE ✓")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Proteomics Pipeline — Step 2: HCP Analysis"
    )
    parser.add_argument("--no-s3", action="store_true",
                        help="Skip S3 upload")
    args = parser.parse_args()

    logger.remove()
    logger.add(sys.stderr, level="INFO", colorize=True)
    run_analysis(upload_s3=not args.no_s3)
