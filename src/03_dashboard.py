"""
src/03_dashboard.py
────────────────────
Step 3 of the Proteomics HCP Pipeline.

Interactive Streamlit dashboard for HCP analysis visualization.
Designed to look like a real biopharma analytical tool.

Run:
  streamlit run src/03_dashboard.py
"""

import json
import sqlite3
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

import streamlit as st

# ── Page config ───────────────────────────────────────────────
st.set_page_config(
    page_title="HCP Analytics Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Paths ─────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).resolve().parent.parent
DATA_PROC   = BASE_DIR / "data/processed"
RESULTS_DIR = BASE_DIR / "results"
DB_PATH     = DATA_PROC / "proteomics.db"

# ── Color palette ─────────────────────────────────────────────
RISK_COLORS = {
    "High":   "#ef4444",
    "Medium": "#f97316",
    "Low":    "#eab308",
    "Trace":  "#6b7280",
}


# ── Data Loading ──────────────────────────────────────────────

@st.cache_data
def load_data():
    proteins_path = DATA_PROC / "proteins_quantified.csv"
    psms_path     = DATA_PROC / "psms_clean.csv"
    report_path   = RESULTS_DIR / "hcp_analysis_report.json"

    if not proteins_path.exists():
        return None, None, None

    proteins = pd.read_csv(proteins_path)
    psms     = pd.read_csv(psms_path)
    report   = json.load(open(report_path)) if report_path.exists() else {}

    return proteins, psms, report


# ── Header ────────────────────────────────────────────────────

st.markdown("""
<style>
    .metric-card {
        background: #1e293b;
        border-radius: 8px;
        padding: 16px;
        border-left: 4px solid #3b82f6;
    }
    .risk-high   { color: #ef4444; font-weight: bold; }
    .risk-medium { color: #f97316; font-weight: bold; }
    .risk-low    { color: #eab308; font-weight: bold; }
</style>
""", unsafe_allow_html=True)

st.title("🧬 Host Cell Protein (HCP) Analytics Dashboard")
st.markdown(
    "**Dataset:** PXD000509 — E. coli HCP Identification in Therapeutic Protein Purification  \n"
    "**Instrument:** LTQ Orbitrap Velos | **Search Engine:** Sequest / Proteome Discoverer v1.3  \n"
    "**Reference:** Bomans et al., Roche Diagnostics GmbH"
)
st.divider()

# ── Load data ─────────────────────────────────────────────────
proteins, psms, report = load_data()

if proteins is None:
    st.error(
        "⚠️ No processed data found. Please run the pipeline first:\n\n"
        "```bash\n"
        "python src/01_parse_mztab.py\n"
        "python src/02_hcp_analysis.py --no-s3\n"
        "```"
    )
    st.stop()

# ── Sidebar filters ───────────────────────────────────────────
with st.sidebar:
    st.header("🔍 Filters")

    min_psms = st.slider(
        "Minimum PSMs",
        min_value=1,
        max_value=int(proteins["total_psms"].max()) if "total_psms" in proteins.columns else 10,
        value=1,
    )

    if "hcp_risk" in proteins.columns:
        risk_filter = st.multiselect(
            "HCP Risk Level",
            options=["High", "Medium", "Low", "Trace"],
            default=["High", "Medium", "Low", "Trace"],
        )
    else:
        risk_filter = []

    st.divider()
    st.markdown("**About**")
    st.markdown(
        "This dashboard analyzes host cell protein (HCP) impurities "
        "in biopharmaceutical purification using LC-MS/MS proteomics data. "
        "Spectral counting is used as a semi-quantitative abundance measure."
    )

# ── Apply filters ─────────────────────────────────────────────
df = proteins.copy()
if "total_psms" in df.columns:
    df = df[df["total_psms"] >= min_psms]
if risk_filter and "hcp_risk" in df.columns:
    df = df[df["hcp_risk"].isin(risk_filter)]

# ── KPI Metrics Row ───────────────────────────────────────────
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    st.metric("HCPs Identified", len(proteins))
with col2:
    st.metric("After Filters", len(df))
with col3:
    total_psms = report.get("total_psms", len(psms))
    st.metric("Total PSMs", f"{total_psms:,}")
with col4:
    unique_pep = report.get("unique_peptides", psms["sequence"].nunique() if "sequence" in psms.columns else 0)
    st.metric("Unique Peptides", f"{unique_pep:,}")
with col5:
    mass_acc = report.get("mass_accuracy", {})
    mean_ppm = mass_acc.get("mean_ppm", None)
    st.metric("Mean Mass Error", f"{mean_ppm:.3f} ppm" if mean_ppm else "N/A")

st.divider()

# ── Main Charts ───────────────────────────────────────────────
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 HCP Abundance", "🎯 Risk Analysis", "🔬 PSM Quality", "📋 Data Table"
])

with tab1:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("Top HCPs by Spectral Count")
        top_n = st.slider("Show top N proteins", 5, min(30, len(df)), 15)
        plot_df = df.nlargest(top_n, "total_psms") if "total_psms" in df.columns else df.head(top_n)

        label_col = "gene_name" if "gene_name" in plot_df.columns else "uniprot_accession"
        plot_df["label"] = plot_df[label_col].fillna(plot_df["uniprot_accession"])

        color_col = "hcp_risk" if "hcp_risk" in plot_df.columns else None
        fig = px.bar(
            plot_df.sort_values("total_psms"),
            x="total_psms",
            y="label",
            orientation="h",
            color=color_col,
            color_discrete_map=RISK_COLORS if color_col else None,
            labels={"total_psms": "PSM Count", "label": "Protein"},
            title=f"Top {top_n} Most Abundant HCPs (Spectral Counting)",
            template="plotly_dark",
        )
        fig.update_layout(height=500, showlegend=True)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Spectral Abundance %")
        if "spectral_abundance_pct" in df.columns:
            top_pie = df.nlargest(8, "total_psms").copy()
            label_col2 = "gene_name" if "gene_name" in top_pie.columns else "uniprot_accession"
            top_pie["label"] = top_pie[label_col2].fillna(top_pie["uniprot_accession"])

            others_pct = 100 - top_pie["spectral_abundance_pct"].sum()
            if others_pct > 0:
                others_row = pd.DataFrame([{"label": "Others", "spectral_abundance_pct": round(others_pct, 2)}])
                top_pie = pd.concat([top_pie, others_row], ignore_index=True)

            fig2 = px.pie(
                top_pie,
                values="spectral_abundance_pct",
                names="label",
                title="Spectral Abundance Distribution",
                template="plotly_dark",
                hole=0.4,
            )
            fig2.update_layout(height=400)
            st.plotly_chart(fig2, use_container_width=True)

            # Key finding callout
            top_hcp = df.nlargest(1, "total_psms").iloc[0]
            top_label = top_hcp.get("gene_name", top_hcp.get("uniprot_accession", ""))
            st.info(
                f"**Key Finding:** {top_label} is the most abundant HCP "
                f"({top_hcp.get('spectral_abundance_pct', 0):.1f}% of total spectral counts). "
                f"Bacterial alkaline phosphatase (phoA) was identified as the dominant contaminant "
                f"across purification steps in this study."
            )

with tab2:
    st.subheader("HCP Risk Assessment")
    col1, col2 = st.columns(2)

    with col1:
        if "hcp_risk" in proteins.columns:
            risk_counts = proteins["hcp_risk"].value_counts().reset_index()
            risk_counts.columns = ["Risk Level", "Count"]
            fig3 = px.bar(
                risk_counts,
                x="Risk Level",
                y="Count",
                color="Risk Level",
                color_discrete_map=RISK_COLORS,
                title="HCP Risk Level Distribution",
                template="plotly_dark",
            )
            st.plotly_chart(fig3, use_container_width=True)

    with col2:
        if "total_psms" in df.columns and "hcp_risk" in df.columns:
            label_col3 = "gene_name" if "gene_name" in df.columns else "uniprot_accession"
            df["label"] = df[label_col3].fillna(df["uniprot_accession"])

            fig4 = px.scatter(
                df,
                x="num_peptides_distinct" if "num_peptides_distinct" in df.columns else "total_psms",
                y="total_psms",
                color="hcp_risk",
                color_discrete_map=RISK_COLORS,
                hover_data=["label", "uniprot_accession"],
                size="total_psms",
                title="PSMs vs Unique Peptides by Risk Level",
                template="plotly_dark",
                labels={
                    "num_peptides_distinct": "Unique Peptides",
                    "total_psms": "PSM Count",
                },
            )
            st.plotly_chart(fig4, use_container_width=True)

    # High risk table
    st.subheader("⚠️ High Risk HCPs")
    high_risk = proteins[proteins["hcp_risk"] == "High"] if "hcp_risk" in proteins.columns else pd.DataFrame()
    if not high_risk.empty:
        display_cols = [c for c in ["uniprot_accession", "gene_name", "protein_name",
                                     "total_psms", "spectral_abundance_pct", "organism"]
                        if c in high_risk.columns]
        st.dataframe(high_risk[display_cols], use_container_width=True)
    else:
        st.info("No high-risk HCPs identified with current filters.")

with tab3:
    st.subheader("PSM Quality Metrics")
    col1, col2 = st.columns(2)

    with col1:
        if "mass_error_ppm" in psms.columns:
            errors = psms["mass_error_ppm"].dropna()
            fig5 = px.histogram(
                psms.dropna(subset=["mass_error_ppm"]),
                x="mass_error_ppm",
                nbins=50,
                title="Mass Error Distribution (ppm)",
                template="plotly_dark",
                color_discrete_sequence=["#3b82f6"],
                labels={"mass_error_ppm": "Mass Error (ppm)"},
            )
            fig5.add_vline(x=0, line_dash="dash", line_color="white", annotation_text="0 ppm")
            fig5.add_vline(x=5, line_dash="dot", line_color="#ef4444", annotation_text="+5 ppm")
            fig5.add_vline(x=-5, line_dash="dot", line_color="#ef4444", annotation_text="-5 ppm")
            st.plotly_chart(fig5, use_container_width=True)

            ma = report.get("mass_accuracy", {})
            if ma:
                st.metric("Mean Error", f"{ma.get('mean_ppm', 0):.3f} ppm")
                st.metric("Within ±5 ppm", f"{ma.get('within_5ppm', 0)}%")

    with col2:
        if "charge" in psms.columns:
            charge_df = psms["charge"].dropna().astype(int).value_counts().reset_index()
            charge_df.columns = ["Charge State", "Count"]
            fig6 = px.bar(
                charge_df.sort_values("Charge State"),
                x="Charge State",
                y="Count",
                title="Charge State Distribution",
                template="plotly_dark",
                color_discrete_sequence=["#8b5cf6"],
            )
            st.plotly_chart(fig6, use_container_width=True)

        if "peptide_length" in psms.columns:
            fig7 = px.histogram(
                psms.dropna(subset=["peptide_length"]),
                x="peptide_length",
                nbins=30,
                title="Peptide Length Distribution",
                template="plotly_dark",
                color_discrete_sequence=["#10b981"],
                labels={"peptide_length": "Peptide Length (aa)"},
            )
            st.plotly_chart(fig7, use_container_width=True)

with tab4:
    st.subheader("📋 Protein Data Table")

    search = st.text_input("Search by gene name, accession, or protein name")
    display_df = df.copy()
    if search:
        mask = pd.Series(False, index=display_df.index)
        for col in ["gene_name", "uniprot_accession", "protein_name"]:
            if col in display_df.columns:
                mask |= display_df[col].astype(str).str.contains(search, case=False, na=False)
        display_df = display_df[mask]

    display_cols = [c for c in [
        "abundance_rank", "uniprot_accession", "gene_name", "protein_name",
        "organism", "total_psms", "unique_peptides", "spectral_abundance_pct",
        "avg_mass_error_ppm", "hcp_risk"
    ] if c in display_df.columns]

    st.dataframe(
        display_df[display_cols].reset_index(drop=True),
        use_container_width=True,
        height=400,
    )

    # Download button
    csv = display_df[display_cols].to_csv(index=False)
    st.download_button(
        "⬇️ Download filtered results as CSV",
        data=csv,
        file_name="hcp_results_filtered.csv",
        mime="text/csv",
    )

# ── Footer ────────────────────────────────────────────────────
st.divider()
st.caption(
    "HCP Analytics Pipeline | PXD000509 | "
    "Built with Python, pyteomics, Streamlit, and Plotly | "
    "Data: PRIDE Archive, EBI"
)
