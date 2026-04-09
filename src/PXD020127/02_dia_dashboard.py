"""
src/PXD020127/02_dia_dashboard.py
───────────────────────────────────
Part 2, Step 2 — DIA HCP Analytics Dashboard

Extends the Part 1 dashboard with:
- Multi-condition abundance comparison
- Fold change volcano plots
- DIA reproducibility (CV) analysis
- Cross-project DDA vs DIA comparison tab

Run:
  streamlit run src/PXD020127/02_dia_dashboard.py
"""

import io
import json
import boto3
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import streamlit as st

st.set_page_config(
    page_title="CHO HCP DIA Dashboard — PXD020127",
    page_icon="🔬",
    layout="wide",
)

BASE_DIR    = Path(__file__).resolve().parent.parent.parent
PROC_DIR    = BASE_DIR / "data/processed/PXD020127"
RESULTS_DIR = BASE_DIR / "results/PXD020127"

S3_BUCKET  = "healthcare-etl-mengqi"
S3_PREFIX  = "proteomics/PXD020127"
AWS_REGION = "us-east-1"

RISK_COLORS = {
    "High":   "#ef4444",
    "Medium": "#f97316",
    "Low":    "#eab308",
    "Trace":  "#6b7280",
}

CONDITION_COLORS = {
    "control":     "#3b82f6",
    "lowDO":       "#ef4444",
    "lowTEMP":     "#8b5cf6",
    "low_TEMP_DO": "#f97316",
}

CONDITION_LABELS = {
    "control":     "Control",
    "lowDO":       "Low DO",
    "lowTEMP":     "Low Temp",
    "low_TEMP_DO": "Low Temp+DO",
}


@st.cache_data(ttl=3600)
def load_data():
    """Load DIA data from S3 or local."""
    try:
        s3 = boto3.client("s3", region_name=AWS_REGION)

        df = pd.read_csv(io.BytesIO(
            s3.get_object(Bucket=S3_BUCKET,
                         Key=f"{S3_PREFIX}/processed/proteins_dia_quantified.csv")["Body"].read()
        ))
        df_long = pd.read_csv(io.BytesIO(
            s3.get_object(Bucket=S3_BUCKET,
                         Key=f"{S3_PREFIX}/processed/proteins_dia_long.csv")["Body"].read()
        ))
        summary = json.loads(
            s3.get_object(Bucket=S3_BUCKET,
                         Key=f"{S3_PREFIX}/results/dia_analysis_summary.json")["Body"].read()
        )
        return df, df_long, summary, "s3"
    except Exception:
        try:
            df      = pd.read_csv(PROC_DIR / "proteins_dia_quantified.csv")
            df_long = pd.read_csv(PROC_DIR / "proteins_dia_long.csv")
            summary = json.load(open(RESULTS_DIR / "dia_analysis_summary.json"))
            return df, df_long, summary, "local"
        except Exception:
            return None, None, None, "none"


# ── Header ────────────────────────────────────────────────────
st.title("🔬 CHO HCP Analytics Dashboard — DIA LC-MS/MS")
st.markdown(
    "**Dataset:** PXD020127 — CHO Cell HCP Profiling Across Bioprocess Conditions  \n"
    "**Acquisition:** Data-Independent Acquisition (DIA) | **Host:** CHO (*Cricetulus griseus*)  \n"
    "**Conditions:** Control · Low DO · Low Temperature · Low Temp+DO  \n"
    "**Reference:** [PRIDE Archive PXD020127](https://www.ebi.ac.uk/pride/archive/projects/PXD020127)"
)
st.divider()

with st.spinner("Loading DIA data..."):
    df, df_long, summary, source = load_data()

if df is None:
    st.error("Run `python src/PXD020127/01_parse_dia.py` first.")
    st.stop()

if source == "s3":
    st.success(f"✅ Live data from AWS S3 — s3://{S3_BUCKET}/{S3_PREFIX}/")
else:
    st.info("📁 Local data")

# ── Sidebar ───────────────────────────────────────────────────
with st.sidebar:
    st.header("🔍 Filters")
    risk_filter = st.multiselect(
        "HCP Risk Level",
        ["High", "Medium", "Low", "Trace"],
        default=["High", "Medium", "Low", "Trace"]
    )
    min_abundance = st.number_input(
        "Min mean abundance", value=0, min_value=0
    )
    show_responsive_only = st.checkbox("Show bioprocess-responsive only", False)
    st.divider()
    st.markdown("**Dataset Info**")
    st.markdown(f"**Proteins:** {len(df):,}")
    st.markdown(f"**Conditions:** 4 bioprocess")
    st.markdown(f"**Replicates:** 3 per condition")
    st.markdown(f"**Source:** `{source.upper()}`")

# Apply filters
df_filt = df.copy()
if risk_filter and "hcp_risk" in df_filt.columns:
    df_filt = df_filt[df_filt["hcp_risk"].isin(risk_filter)]
if "overall_mean_abundance" in df_filt.columns:
    df_filt = df_filt[df_filt["overall_mean_abundance"] >= min_abundance]
if show_responsive_only and "num_conditions_responsive" in df_filt.columns:
    df_filt = df_filt[df_filt["num_conditions_responsive"] >= 1]

# ── KPI Row ───────────────────────────────────────────────────
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Total HCPs", f"{len(df):,}")
c2.metric("After Filters", f"{len(df_filt):,}")
c3.metric("High Risk", len(df[df["hcp_risk"] == "High"]) if "hcp_risk" in df.columns else "N/A")
c4.metric("Stress-Responsive", int(df["num_conditions_responsive"].ge(2).sum()) if "num_conditions_responsive" in df.columns else "N/A")
c5.metric("Conditions", "4 bioprocess")
st.divider()

# ── Tabs ──────────────────────────────────────────────────────
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Condition Comparison",
    "🌋 Fold Change Analysis",
    "📈 Reproducibility (CV)",
    "⚖️ DDA vs DIA",
    "📋 Data Table",
])

with tab1:
    st.subheader("HCP Abundance Across Bioprocess Conditions")
    col1, col2 = st.columns([2, 1])

    with col1:
        top_n = st.slider("Top N proteins", 10, min(50, len(df_filt)), 20, key="cond_topn")
        top_df = df_filt.nsmallest(top_n, "abundance_rank").copy()

        mean_cols = {cond: f"{cond}_mean" for cond in CONDITION_LABELS
                    if f"{cond}_mean" in top_df.columns}
        label_col = "gene_name" if "gene_name" in top_df.columns else "accession"
        top_df["label"] = top_df[label_col].fillna(top_df["accession"])

        # Melt for grouped bar chart
        melt_df = top_df[["label"] + list(mean_cols.values())].melt(
            id_vars="label",
            value_vars=list(mean_cols.values()),
            var_name="condition_raw",
            value_name="mean_abundance",
        )
        melt_df["condition"] = melt_df["condition_raw"].str.replace("_mean", "")
        melt_df["condition_label"] = melt_df["condition"].map(CONDITION_LABELS)

        fig = px.bar(
            melt_df,
            x="label", y="mean_abundance",
            color="condition_label",
            barmode="group",
            title=f"Top {top_n} HCPs — Abundance by Bioprocess Condition",
            template="plotly_dark",
            labels={"mean_abundance": "Mean Abundance", "label": "Protein"},
            color_discrete_sequence=list(CONDITION_COLORS.values()),
        )
        fig.update_layout(height=500, xaxis_tickangle=-45)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Risk Distribution")
        if "hcp_risk" in df.columns:
            risk_df = df["hcp_risk"].value_counts().reset_index()
            risk_df.columns = ["Risk", "Count"]
            fig2 = px.pie(
                risk_df, values="Count", names="Risk",
                color="Risk", color_discrete_map=RISK_COLORS,
                hole=0.4, template="plotly_dark",
                title="HCP Risk Distribution",
            )
            st.plotly_chart(fig2, use_container_width=True)

        # Key finding
        if "hcp_risk" in df.columns:
            high_risk = df[df["hcp_risk"] == "High"]
            st.warning(
                f"**{len(high_risk)} High-Risk HCPs identified** across 4 bioprocess conditions. "
                f"Stress conditions (low DO, low temperature) alter HCP profiles, "
                f"potentially introducing risk-relevant proteins not detected under control conditions."
            )

with tab2:
    st.subheader("🌋 Fold Change vs Control — Volcano Plots")
    stress_cond = st.selectbox(
        "Select stress condition",
        ["lowDO", "lowTEMP", "low_TEMP_DO"],
        format_func=lambda x: CONDITION_LABELS[x],
    )

    fc_col = f"{stress_cond}_log2fc"
    if fc_col in df_filt.columns:
        plot_df = df_filt.dropna(subset=[fc_col, "overall_mean_abundance"]).copy()
        plot_df["label"] = plot_df.get("gene_name", plot_df["accession"]).fillna(plot_df["accession"])
        plot_df["-log10_abundance"] = np.log10(plot_df["overall_mean_abundance"] + 1)
        plot_df["significant"] = plot_df[fc_col].abs() > 1

        fig3 = px.scatter(
            plot_df,
            x=fc_col,
            y="-log10_abundance",
            color="hcp_risk",
            color_discrete_map=RISK_COLORS,
            hover_data=["label", "accession", "unique_peptides"],
            title=f"Fold Change: {CONDITION_LABELS[stress_cond]} vs Control",
            template="plotly_dark",
            labels={
                fc_col: "Log2 Fold Change vs Control",
                "-log10_abundance": "Log10 Mean Abundance",
            },
        )
        fig3.add_vline(x=1,  line_dash="dash", line_color="#ef4444", annotation_text="+2x")
        fig3.add_vline(x=-1, line_dash="dash", line_color="#3b82f6", annotation_text="-2x")
        fig3.add_vline(x=0,  line_dash="dot",  line_color="white")
        fig3.update_layout(height=500)
        st.plotly_chart(fig3, use_container_width=True)

        up   = (plot_df[fc_col] > 1).sum()
        down = (plot_df[fc_col] < -1).sum()
        c1, c2, c3 = st.columns(3)
        c1.metric("Upregulated (>2x)", up)
        c2.metric("Downregulated (<0.5x)", down)
        c3.metric("Unchanged", len(plot_df) - up - down)

with tab3:
    st.subheader("📈 DIA Reproducibility — Coefficient of Variation")
    st.markdown("CV < 20% is the industry benchmark for acceptable DIA quantification reproducibility.")

    cv_data = []
    for cond, label in CONDITION_LABELS.items():
        cv_col = f"{cond}_cv_pct"
        if cv_col in df.columns:
            vals = df[cv_col].dropna()
            cv_data.append({
                "Condition":     label,
                "Median CV (%)": round(vals.median(), 2),
                "Mean CV (%)":   round(vals.mean(), 2),
                "% < 20%":       round((vals < 20).mean() * 100, 1),
                "% < 10%":       round((vals < 10).mean() * 100, 1),
            })

    if cv_data:
        cv_summary = pd.DataFrame(cv_data)
        st.dataframe(cv_summary, use_container_width=True)

        # CV distribution plot
        cv_cols = {cond: f"{cond}_cv_pct" for cond in CONDITION_LABELS
                  if f"{cond}_cv_pct" in df.columns}
        cv_melt = df[list(cv_cols.values())].melt(
            var_name="condition_raw", value_name="cv_pct"
        ).dropna()
        cv_melt["condition"] = cv_melt["condition_raw"].str.replace("_cv_pct", "")
        cv_melt["condition_label"] = cv_melt["condition"].map(CONDITION_LABELS)

        fig5 = px.box(
            cv_melt, x="condition_label", y="cv_pct",
            color="condition_label",
            color_discrete_sequence=list(CONDITION_COLORS.values()),
            title="CV Distribution by Condition",
            template="plotly_dark",
            labels={"cv_pct": "CV (%)", "condition_label": "Condition"},
        )
        fig5.add_hline(y=20, line_dash="dash", line_color="#ef4444",
                       annotation_text="20% threshold")
        fig5.update_layout(height=400, showlegend=False)
        st.plotly_chart(fig5, use_container_width=True)

with tab4:
    st.subheader("⚖️ DDA vs DIA — Technology Comparison")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("### PXD000509 — DDA (E. coli)")
        st.markdown("""
        - **Mode:** Data-Dependent Acquisition
        - **Host:** *Escherichia coli*
        - **Proteins:** 54 HCPs identified
        - **PSMs:** 220 total
        - **Quantification:** Spectral counting (semi-quantitative)
        - **Key HCP:** phoA (Bacterial alkaline phosphatase)
        - **Mass accuracy:** -1.24 ppm (LTQ Orbitrap Velos)
        - **Limitation:** Stochastic sampling — misses low-abundance HCPs
        """)

    with col2:
        st.markdown("### PXD020127 — DIA (CHO)")
        st.markdown(f"""
        - **Mode:** Data-Independent Acquisition
        - **Host:** *Cricetulus griseus* (CHO)
        - **Proteins:** {len(df):,} HCPs quantified
        - **Conditions:** 4 bioprocess conditions × 3 replicates
        - **Quantification:** Area-under-curve (quantitative)
        - **Key advantage:** Comprehensive, reproducible, condition-comparable
        - **Bioprocess insight:** Stress conditions alter HCP profiles
        - **Industry trend:** Emerging gold standard for HCP monitoring
        """)

    st.divider()
    comparison_data = {
        "Feature":            ["Coverage", "Reproducibility", "Quantification", "Throughput", "Sensitivity"],
        "DDA (PXD000509)":    ["Stochastic", "Variable", "Spectral counting", "High", "Good"],
        "DIA (PXD020127)":    ["Comprehensive", "High (CV<20%)", "Area-under-curve", "Medium", "Excellent"],
    }
    st.dataframe(pd.DataFrame(comparison_data), use_container_width=True, hide_index=True)

    st.info(
        "**Portfolio insight:** Demonstrating proficiency in both DDA and DIA workflows "
        "covers the full spectrum of current LC-MS/MS HCP analysis approaches used in "
        "biopharmaceutical development. DDA is the established method; DIA is the emerging "
        "standard for comprehensive, quantitative HCP monitoring."
    )

with tab5:
    st.subheader("📋 Full Protein Data Table")
    search = st.text_input("Search by gene name or accession")
    disp = df_filt.copy()
    if search:
        mask = pd.Series(False, index=disp.index)
        for col in ["gene_name", "accession", "description"]:
            if col in disp.columns:
                mask |= disp[col].astype(str).str.contains(search, case=False, na=False)
        disp = disp[mask]

    show_cols = [c for c in [
        "abundance_rank", "accession", "gene_name", "unique_peptides",
        "overall_mean_abundance", "control_mean", "lowDO_mean",
        "lowTEMP_mean", "low_TEMP_DO_mean", "lowDO_log2fc",
        "lowTEMP_log2fc", "num_conditions_responsive", "hcp_risk",
    ] if c in disp.columns]

    st.dataframe(disp[show_cols].reset_index(drop=True), use_container_width=True, height=400)
    st.download_button(
        "⬇️ Download CSV",
        data=disp[show_cols].to_csv(index=False),
        file_name="cho_hcp_dia_results.csv",
        mime="text/csv",
    )

st.divider()
c1, c2, c3 = st.columns(3)
c1.caption("🔬 CHO HCP DIA Pipeline | PXD020127")
c2.caption(f"☁️ AWS S3: {S3_BUCKET}")
c3.caption("Python · Streamlit · Plotly · pandas")
