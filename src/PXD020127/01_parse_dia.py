"""
src/PXD020127/01_parse_dia.py
──────────────────────────────
Part 2, Step 1 of the Proteomics HCP Pipeline.

Dataset: PXD020127 — CHO cell HCP analysis using DIA LC-MS/MS
         Automated magnetic bead-based sample prep + DIA acquisition
         4 bioprocess conditions x 3 replicates = 12 samples

What this script does:
  1. Loads Results.xlsx (1,184 proteins x 30 columns)
  2. Cleans and restructures into long format
  3. Calculates per-condition statistics (mean, CV)
  4. Classifies HCP risk based on known CHO HCPs
  5. Compares conditions to identify bioprocess-responsive HCPs
  6. Saves to CSV and SQLite

Bioprocess conditions:
  control     — standard CHO culture conditions
  lowDO       — low dissolved oxygen stress
  lowTEMP     — low temperature (hypothermic) culture
  low_TEMP_DO — combined low temperature + low DO stress

Run:
  python src/PXD020127/01_parse_dia.py
"""

import sys
import json
import sqlite3
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timezone
from loguru import logger

# ── Paths ─────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).resolve().parent.parent.parent
RAW_DIR     = BASE_DIR / "data/raw/PXD020127"
PROC_DIR    = BASE_DIR / "data/processed/PXD020127"
RESULTS_DIR = BASE_DIR / "results/PXD020127"
DB_PATH     = BASE_DIR / "data/processed/proteomics.db"

PROC_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

EXCEL_FILE = RAW_DIR / "Results.xlsx"

# ── Sample / Condition mapping ────────────────────────────────
CONDITIONS = {
    "control":     ["control_1",     "control_2",     "control_3"],
    "lowDO":       ["lowDO_1",       "lowDO_2",       "lowDO_3"],
    "lowTEMP":     ["lowTEMP_1",     "lowTEMP_2",     "lowTEMP_3"],
    "low_TEMP_DO": ["low_TEMP_DO_1", "low_TEMP_DO_2", "low_TEMP_DO_3"],
}

PPM_COLS = {
    "control":     ["control_1 ppm",     "control_2 ppm",     "control_3 ppm"],
    "lowDO":       ["lowDO_1 ppm",       "lowDO_2 ppm",       "lowDO_3 ppm"],
    "lowTEMP":     ["lowTEMP_1 ppm",     "lowTEMP_2 ppm",     "lowTEMP_3 ppm"],
    "low_TEMP_DO": ["low_TEMP_DO_1 ppm", "low_TEMP_DO_2 ppm", "low_TEMP_DO_3 ppm"],
}

CONDITION_LABELS = {
    "control":     "Control (Standard)",
    "lowDO":       "Low Dissolved Oxygen",
    "lowTEMP":     "Low Temperature",
    "low_TEMP_DO": "Low Temp + Low DO",
}

# Known high-risk CHO HCPs from biopharma literature
HIGH_RISK_CHO_GENES = {
    "plbl2",   # Phospholipase B-like 2 — most studied CHO HCP
    "lpl",     # Lipoprotein lipase — lipase activity risk
    "pla2g15", # Phospholipase A2 — lipase activity risk
    "cathb",   # Cathepsin B — protease risk
    "ctsb",    # Cathepsin B alternate name
    "ctsd",    # Cathepsin D — protease
    "hcpcs",   # HCP-related
    "anxa1",   # Annexin A1
    "lgmn",    # Legumain — asparagine endopeptidase
    "eno1",    # Enolase — common contaminant
    "gapdh",   # GAPDH — common contaminant
    "hsp90aa1",# Heat shock protein
    "hspa5",   # BiP/GRP78 — ER chaperone
    "hspa8",   # Heat shock cognate 70
}


# ── Data Loading & Cleaning ───────────────────────────────────

def load_and_clean(filepath: Path) -> pd.DataFrame:
    """Load Excel results and clean column names."""
    logger.info(f"Loading: {filepath.name}")
    df = pd.read_excel(filepath, sheet_name="Tabelle1")
    logger.info(f"Loaded {len(df):,} proteins x {len(df.columns)} columns")

    # Clean column names
    df = df.rename(columns={
        "Accession":       "accession",
        "Description":     "description",
        "Genenames":       "gene_name",
        "Unique peptides": "unique_peptides",
        "Theoretical Mass ANTIBODY (Dalton)": "antibody_mass_da",
        "MASS (Dalton)":   "protein_mass_da",
    })

    # Cast abundance columns to numeric
    all_sample_cols = [c for cond in CONDITIONS.values() for c in cond]
    all_ppm_cols    = [c for cols in PPM_COLS.values() for c in cols]

    for col in all_sample_cols + all_ppm_cols + ["unique_peptides", "protein_mass_da"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Extract species from description
    df["species"] = df["description"].str.extract(r"OS=([^O]+?)(?:\s+OX=|$)")
    df["species"]  = df["species"].str.strip()

    # Clean gene names
    df["gene_name"] = df["gene_name"].astype(str).str.strip()
    df["gene_name_lower"] = df["gene_name"].str.lower()

    logger.info(f"Species found: {df['species'].value_counts().head(3).to_dict()}")
    return df


# ── Per-Condition Statistics ──────────────────────────────────

def calculate_condition_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate mean abundance, CV (coefficient of variation),
    and detection rate per bioprocess condition.
    CV = std/mean * 100 — key DIA quality metric, should be < 20%
    """
    logger.info("Calculating per-condition statistics...")

    for cond, cols in CONDITIONS.items():
        available = [c for c in cols if c in df.columns]
        if not available:
            continue

        df[f"{cond}_mean"]     = df[available].mean(axis=1)
        df[f"{cond}_std"]      = df[available].std(axis=1)
        df[f"{cond}_cv_pct"]   = (df[f"{cond}_std"] / df[f"{cond}_mean"] * 100).round(2)
        df[f"{cond}_detected"] = (df[available] > 0).sum(axis=1)

        # Mean ppm mass accuracy per condition
        ppm_cols = [c for c in PPM_COLS.get(cond, []) if c in df.columns]
        if ppm_cols:
            df[f"{cond}_mean_ppm"] = df[ppm_cols].mean(axis=1).round(4)

    # Overall mean across all conditions
    all_means = [f"{c}_mean" for c in CONDITIONS if f"{c}_mean" in df.columns]
    df["overall_mean_abundance"] = df[all_means].mean(axis=1)
    df["abundance_rank"] = df["overall_mean_abundance"].rank(ascending=False, method="min").astype(int)

    return df


# ── Fold Change Analysis ──────────────────────────────────────

def calculate_fold_changes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate fold changes vs control for each stress condition.
    This identifies bioprocess-responsive HCPs — proteins whose
    abundance changes under different culture conditions.
    log2 fold change > 1 = 2x upregulation
    log2 fold change < -1 = 2x downregulation
    """
    logger.info("Calculating fold changes vs control...")

    stress_conditions = ["lowDO", "lowTEMP", "low_TEMP_DO"]

    for cond in stress_conditions:
        if f"{cond}_mean" not in df.columns or "control_mean" not in df.columns:
            continue

        # Log2 fold change
        df[f"{cond}_log2fc"] = np.log2(
            (df[f"{cond}_mean"] + 1) / (df["control_mean"] + 1)
        ).round(4)

        # Significance flag (simple threshold — |log2FC| > 1)
        df[f"{cond}_responsive"] = df[f"{cond}_log2fc"].abs() > 1

    # Count how many conditions each protein responds to
    responsive_cols = [f"{c}_responsive" for c in stress_conditions
                      if f"{c}_responsive" in df.columns]
    if responsive_cols:
        df["num_conditions_responsive"] = df[responsive_cols].sum(axis=1)

    return df


# ── HCP Risk Classification ───────────────────────────────────

def classify_hcp_risk(df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify CHO HCP risk level based on:
    1. Known high-risk gene names from literature
    2. Overall abundance rank
    3. Bioprocess responsiveness (stress-induced proteins are higher risk)
    """
    def get_risk(row):
        gene  = str(row.get("gene_name_lower", "")).lower()
        rank  = row.get("abundance_rank", 9999)
        responsive = row.get("num_conditions_responsive", 0)

        if gene in HIGH_RISK_CHO_GENES:
            return "High"
        elif rank <= 50 and responsive >= 2:
            return "High"
        elif rank <= 100 or responsive >= 2:
            return "Medium"
        elif rank <= 300:
            return "Low"
        else:
            return "Trace"

    df["hcp_risk"] = df.apply(get_risk, axis=1)
    return df


# ── Comparison with PXD000509 ─────────────────────────────────

def cross_project_summary(df: pd.DataFrame) -> dict:
    """
    Generate a comparison summary between this DIA dataset
    and the DDA dataset from Part 1 (PXD000509).
    """
    return {
        "project":           "PXD020127",
        "acquisition_mode":  "DIA (Data-Independent Acquisition)",
        "host_cell":         "CHO (Chinese Hamster Ovary)",
        "total_proteins":    len(df),
        "conditions":        list(CONDITION_LABELS.values()),
        "replicates_per_condition": 3,
        "comparison_with_pxd000509": {
            "PXD000509": {
                "mode":       "DDA",
                "host":       "E. coli",
                "proteins":   54,
                "psms":       220,
            },
            "PXD020127": {
                "mode":       "DIA",
                "host":       "CHO",
                "proteins":   len(df),
                "conditions": 4,
            },
            "key_difference": (
                "DIA provides comprehensive, unbiased quantification across all "
                "conditions simultaneously. DDA stochastically samples the most "
                "abundant precursors. DIA is emerging as the gold standard for "
                "HCP monitoring in biopharmaceutical QC workflows."
            ),
        },
    }


# ── Storage ───────────────────────────────────────────────────

def save_results(df: pd.DataFrame, summary: dict):
    """Save processed data to CSV and SQLite."""
    # Save main protein table
    out_path = PROC_DIR / "proteins_dia_quantified.csv"
    df.to_csv(out_path, index=False)
    logger.success(f"Saved: {out_path} ({len(df):,} proteins)")

    # Save condition stats in long format for dashboard
    long_rows = []
    for cond, label in CONDITION_LABELS.items():
        if f"{cond}_mean" not in df.columns:
            continue
        for _, row in df.iterrows():
            long_rows.append({
                "accession":    row["accession"],
                "gene_name":    row["gene_name"],
                "condition":    cond,
                "condition_label": label,
                "mean_abundance": row.get(f"{cond}_mean", np.nan),
                "cv_pct":       row.get(f"{cond}_cv_pct", np.nan),
                "log2fc":       row.get(f"{cond}_log2fc", np.nan),
                "detected":     row.get(f"{cond}_detected", np.nan),
                "hcp_risk":     row.get("hcp_risk", "Trace"),
                "abundance_rank": row.get("abundance_rank", 9999),
            })

    df_long = pd.DataFrame(long_rows)
    long_path = PROC_DIR / "proteins_dia_long.csv"
    df_long.to_csv(long_path, index=False)
    logger.success(f"Saved long format: {long_path} ({len(df_long):,} rows)")

    # Save to SQLite
    conn = sqlite3.connect(DB_PATH)
    df.to_sql("proteins_dia", conn, if_exists="replace", index=False)
    df_long.to_sql("proteins_dia_long", conn, if_exists="replace", index=False)
    conn.close()
    logger.success(f"Saved to SQLite: {DB_PATH}")

    # Save summary JSON
    summary_path = RESULTS_DIR / "dia_analysis_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    logger.success(f"Saved summary: {summary_path}")


# ── Print Summary ─────────────────────────────────────────────

def print_summary(df: pd.DataFrame):
    print(f"\n{'='*60}")
    print(f"  DIA HCP ANALYSIS SUMMARY — PXD020127")
    print(f"{'='*60}")
    print(f"  Host cell      : CHO (Cricetulus griseus)")
    print(f"  Acquisition    : DIA (Data-Independent Acquisition)")
    print(f"  Total proteins : {len(df):,}")
    print(f"  Conditions     : {', '.join(CONDITION_LABELS.values())}")

    print(f"\n── Top 10 HCPs by Overall Abundance ──────────────────")
    display_cols = [c for c in ["accession", "gene_name", "overall_mean_abundance",
                                "control_mean", "lowDO_mean", "lowTEMP_mean",
                                "num_conditions_responsive", "hcp_risk"]
                   if c in df.columns]
    print(df.nsmallest(10, "abundance_rank")[display_cols].to_string(index=False))

    print(f"\n── HCP Risk Distribution ──────────────────────────────")
    if "hcp_risk" in df.columns:
        for risk, count in df["hcp_risk"].value_counts().items():
            print(f"  {risk:<10}: {count:,} proteins")

    print(f"\n── Bioprocess-Responsive HCPs ─────────────────────────")
    if "num_conditions_responsive" in df.columns:
        responsive = df[df["num_conditions_responsive"] >= 2]
        print(f"  Proteins responsive in ≥2 stress conditions: {len(responsive)}")
        if len(responsive) > 0:
            resp_cols = [c for c in ["accession", "gene_name", "hcp_risk",
                                     "lowDO_log2fc", "lowTEMP_log2fc",
                                     "num_conditions_responsive"]
                        if c in responsive.columns]
            print(responsive.nlargest(5, "overall_mean_abundance")[resp_cols].to_string(index=False))

    print(f"\n── CV Quality Check (DIA Reproducibility) ─────────────")
    for cond in CONDITIONS:
        cv_col = f"{cond}_cv_pct"
        if cv_col in df.columns:
            median_cv = df[cv_col].median()
            pct_under20 = (df[cv_col] < 20).mean() * 100
            print(f"  {CONDITION_LABELS[cond]:<30}: median CV={median_cv:.1f}%, "
                  f"{pct_under20:.1f}% < 20%")
    print(f"{'='*60}\n")


# ── Main ──────────────────────────────────────────────────────

def run_parse():
    logger.info("=" * 60)
    logger.info("  PART 2 STEP 1 — PARSE DIA RESULTS (PXD020127)")
    logger.info("=" * 60)

    df = load_and_clean(EXCEL_FILE)
    df = calculate_condition_stats(df)
    df = calculate_fold_changes(df)
    df = classify_hcp_risk(df)

    summary = cross_project_summary(df)
    save_results(df, summary)
    print_summary(df)

    logger.success("PART 2 STEP 1 COMPLETE ✓")
    return df, summary


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO", colorize=True)
    run_parse()
