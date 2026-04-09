"""
src/03_dashboard.py - Updated to load from AWS S3 with local fallback
"""
import io
import json
import boto3
import numpy as np
import pandas as pd
import plotly.express as px
from pathlib import Path
import streamlit as st

st.set_page_config(page_title="HCP Analytics Dashboard", page_icon="🧬", layout="wide")

BASE_DIR   = Path(__file__).resolve().parent.parent
DATA_PROC  = BASE_DIR / "data/processed"
RESULTS    = BASE_DIR / "results"
S3_BUCKET  = "healthcare-etl-mengqi"
S3_PREFIX  = "proteomics/PXD000509"
AWS_REGION = "us-east-1"

RISK_COLORS = {"High": "#ef4444", "Medium": "#f97316", "Low": "#eab308", "Trace": "#6b7280"}

@st.cache_data(ttl=3600)
def load_data():
    try:
        s3 = boto3.client("s3", region_name=AWS_REGION)
        proteins = pd.read_csv(io.BytesIO(s3.get_object(Bucket=S3_BUCKET, Key=f"{S3_PREFIX}/processed/proteins_quantified.csv")["Body"].read()))
        psms = pd.read_csv(io.BytesIO(s3.get_object(Bucket=S3_BUCKET, Key=f"{S3_PREFIX}/processed/psms_clean.csv")["Body"].read()))
        report = json.loads(s3.get_object(Bucket=S3_BUCKET, Key=f"{S3_PREFIX}/results/hcp_analysis_report.json")["Body"].read().decode("utf-8"))
        return proteins, psms, report, "s3"
    except Exception:
        try:
            proteins = pd.read_csv(DATA_PROC / "proteins_quantified.csv")
            psms = pd.read_csv(DATA_PROC / "psms_clean.csv")
            report = json.load(open(RESULTS / "hcp_analysis_report.json"))
            return proteins, psms, report, "local"
        except Exception:
            return None, None, None, "none"

st.title("🧬 Host Cell Protein (HCP) Analytics Dashboard")
st.markdown("**Dataset:** PXD000509 — E. coli HCP Identification in Therapeutic Protein Purification  \n**Instrument:** LTQ Orbitrap Velos | **Search Engine:** Sequest / Proteome Discoverer v1.3  \n**Reference:** Bomans et al., Roche Diagnostics GmbH | **Data:** [PRIDE Archive PXD000509](https://www.ebi.ac.uk/pride/archive/projects/PXD000509)")
st.divider()

with st.spinner("Loading from AWS S3..."):
    proteins, psms, report, source = load_data()

if proteins is None:
    st.error("No data found. Run the pipeline first.")
    st.stop()

if source == "s3":
    st.success(f"✅ Live data from AWS S3 — s3://{S3_BUCKET}/{S3_PREFIX}/")
else:
    st.info("📁 Local data")

with st.sidebar:
    st.header("🔍 Filters")
    min_psms = st.slider("Minimum PSMs", 1, int(proteins["total_psms"].max()), 1)
    risk_filter = st.multiselect("HCP Risk Level", ["High","Medium","Low","Trace"], default=["High","Medium","Low","Trace"])
    st.divider()
    st.markdown(f"**Source:** `{source.upper()}`")
    st.markdown(f"**Bucket:** `{S3_BUCKET}`")

df = proteins[proteins["total_psms"] >= min_psms]
if risk_filter:
    df = df[df["hcp_risk"].isin(risk_filter)]

c1,c2,c3,c4,c5 = st.columns(5)
c1.metric("HCPs Identified", len(proteins))
c2.metric("After Filters", len(df))
c3.metric("Total PSMs", f"{report.get('total_psms', len(psms)):,}")
c4.metric("Unique Peptides", f"{report.get('unique_peptides', 0):,}")
mean_ppm = report.get("mass_accuracy", {}).get("mean_ppm")
c5.metric("Mean Mass Error", f"{mean_ppm:.3f} ppm" if mean_ppm else "N/A")
st.divider()

tab1,tab2,tab3,tab4 = st.tabs(["📊 HCP Abundance","🎯 Risk Analysis","🔬 PSM Quality","📋 Data Table"])

with tab1:
    col1,col2 = st.columns([2,1])
    with col1:
        st.subheader("Top HCPs by Spectral Count")
        top_n = st.slider("Show top N", 5, min(30,len(df)), 15)
        plot_df = df.nlargest(top_n,"total_psms").copy()
        plot_df["label"] = plot_df.get("gene_name", plot_df["uniprot_accession"]).fillna(plot_df["uniprot_accession"])
        fig = px.bar(plot_df.sort_values("total_psms"), x="total_psms", y="label", orientation="h",
                     color="hcp_risk", color_discrete_map=RISK_COLORS, template="plotly_dark",
                     labels={"total_psms":"PSM Count","label":"Protein"},
                     title=f"Top {top_n} Most Abundant HCPs (Spectral Counting)")
        fig.update_layout(height=500)
        st.plotly_chart(fig, use_container_width=True)
    with col2:
        st.subheader("Spectral Abundance %")
        top_pie = df.nlargest(8,"total_psms").copy()
        top_pie["label"] = top_pie.get("gene_name", top_pie["uniprot_accession"]).fillna(top_pie["uniprot_accession"])
        others = max(0, 100 - top_pie["spectral_abundance_pct"].sum())
        if others > 0:
            top_pie = pd.concat([top_pie, pd.DataFrame([{"label":"Others","spectral_abundance_pct":round(others,2)}])], ignore_index=True)
        fig2 = px.pie(top_pie, values="spectral_abundance_pct", names="label", hole=0.4, template="plotly_dark", title="Spectral Abundance Distribution")
        fig2.update_layout(height=380)
        st.plotly_chart(fig2, use_container_width=True)
        st.info("**Key Finding:** phoA (P00634, bacterial alkaline phosphatase) is the highest-risk HCP identified — consistent with Bomans et al. Enzymatic BAP activity was monitored via cobas assay across purification steps.")

with tab2:
    col1,col2 = st.columns(2)
    with col1:
        risk_counts = proteins["hcp_risk"].value_counts().reset_index()
        risk_counts.columns = ["Risk Level","Count"]
        fig3 = px.bar(risk_counts, x="Risk Level", y="Count", color="Risk Level", color_discrete_map=RISK_COLORS, template="plotly_dark", title="HCP Risk Distribution")
        st.plotly_chart(fig3, use_container_width=True)
    with col2:
        df2 = df.copy()
        df2["label"] = df2.get("gene_name", df2["uniprot_accession"]).fillna(df2["uniprot_accession"])
        fig4 = px.scatter(df2, x="num_peptides_distinct" if "num_peptides_distinct" in df2.columns else "total_psms",
                          y="total_psms", color="hcp_risk", color_discrete_map=RISK_COLORS,
                          hover_data=["label","uniprot_accession"], size="total_psms",
                          title="PSMs vs Unique Peptides", template="plotly_dark")
        st.plotly_chart(fig4, use_container_width=True)
    st.subheader("⚠️ High Risk HCPs")
    high_risk = proteins[proteins["hcp_risk"]=="High"]
    if not high_risk.empty:
        cols = [c for c in ["uniprot_accession","gene_name","protein_name","total_psms","spectral_abundance_pct","organism"] if c in high_risk.columns]
        st.dataframe(high_risk[cols], use_container_width=True)

with tab3:
    col1,col2 = st.columns(2)
    with col1:
        if "mass_error_ppm" in psms.columns:
            fig5 = px.histogram(psms.dropna(subset=["mass_error_ppm"]), x="mass_error_ppm", nbins=50,
                                title="Mass Error Distribution (ppm)", template="plotly_dark",
                                color_discrete_sequence=["#3b82f6"])
            fig5.add_vline(x=0, line_dash="dash", line_color="white")
            fig5.add_vline(x=5, line_dash="dot", line_color="#ef4444")
            fig5.add_vline(x=-5, line_dash="dot", line_color="#ef4444")
            st.plotly_chart(fig5, use_container_width=True)
            ma = report.get("mass_accuracy", {})
            c1,c2,c3 = st.columns(3)
            c1.metric("Mean Error", f"{ma.get('mean_ppm',0):.3f} ppm")
            c2.metric("Within ±5 ppm", f"{ma.get('within_5ppm',0)}%")
            c3.metric("Std Dev", f"{ma.get('std_ppm',0):.3f} ppm")
    with col2:
        if "charge" in psms.columns:
            charge_df = psms["charge"].dropna().astype(int).value_counts().reset_index()
            charge_df.columns = ["Charge State","Count"]
            fig6 = px.bar(charge_df.sort_values("Charge State"), x="Charge State", y="Count",
                          title="Charge State Distribution", template="plotly_dark",
                          color_discrete_sequence=["#8b5cf6"])
            st.plotly_chart(fig6, use_container_width=True)
        if "peptide_length" in psms.columns:
            fig7 = px.histogram(psms.dropna(subset=["peptide_length"]), x="peptide_length", nbins=30,
                                title="Peptide Length Distribution", template="plotly_dark",
                                color_discrete_sequence=["#10b981"])
            st.plotly_chart(fig7, use_container_width=True)

with tab4:
    st.subheader("📋 Full Protein Data Table")
    search = st.text_input("Search by gene name, accession, or protein name")
    display_df = df.copy()
    if search:
        mask = pd.Series(False, index=display_df.index)
        for col in ["gene_name","uniprot_accession","protein_name"]:
            if col in display_df.columns:
                mask |= display_df[col].astype(str).str.contains(search, case=False, na=False)
        display_df = display_df[mask]
    dcols = [c for c in ["abundance_rank","uniprot_accession","gene_name","protein_name","organism","total_psms","unique_peptides","spectral_abundance_pct","avg_mass_error_ppm","hcp_risk"] if c in display_df.columns]
    st.dataframe(display_df[dcols].reset_index(drop=True), use_container_width=True, height=400)
    st.download_button("⬇️ Download CSV", data=display_df[dcols].to_csv(index=False), file_name="hcp_results.csv", mime="text/csv")

st.divider()
c1,c2,c3 = st.columns(3)
c1.caption("🧬 HCP Analytics Pipeline | PXD000509")
c2.caption(f"☁️ AWS S3: {S3_BUCKET}")
c3.caption("Python · Streamlit · Plotly · pyteomics")
