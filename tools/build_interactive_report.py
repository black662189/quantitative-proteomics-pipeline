#!/usr/bin/env python3
"""Build the offline interactive report from completed analysis tables.
No statistical model is rerun. Plotly is bundled, not loaded from a CDN.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools.report_data import enrich_payload
def records(df):
    return json.loads(df.to_json(orient="records",force_ascii=False))
def build(report_dir,output,plotly_js=None,resource_root=None):
    plotly_js=plotly_js or ROOT/"vendor/plotly-3.1.0.min.js"
    report_dir=Path(report_dir);r=report_dir/"Report"
    cfg=json.loads((r/"7.Provenance/config.resolved.json").read_text(encoding="utf-8"))
    meta=pd.read_csv(r/"1.QC/samples.tsv",sep="\t",keep_default_na=False)
    meta=meta[meta.role=="sample"]
    sample_ids=meta["sample"].tolist()
    ann=pd.read_csv(r/"1.QC/annotation.tsv",sep="\t",keep_default_na=False)
    x=pd.read_csv(r/"5.Expression/normalized_log2.tsv",sep="\t",index_col=0)[sample_ids]
    qc=pd.read_csv(r/"5.Expression/sample_qc.tsv",sep="\t",keep_default_na=False)
    qc=qc[qc.role=="sample"].drop(columns=["column"])
    pcafile=r/"2.Cluster/PCA_scores.tsv"
    pca=pd.read_csv(pcafile,sep="\t") if pcafile.exists() else pd.DataFrame(columns=["sample","PC1","PC2"])
    if len(pca):assert set(pca["sample"])==set(sample_ids)
    full=x.dropna();full=full[full.var(axis=1)>0]
    centered=full.to_numpy().T;centered-=centered.mean(axis=0)
    s=np.linalg.svd(centered,compute_uv=False);explained=(s*s/(s*s).sum())[:2].tolist() if len(s)>=2 and sum(s*s)>0 else []
    payload=dict(project=cfg["project"],species=cfg["species"],taxid=cfg["taxid"],
                 alpha=cfg["alpha"],fc=cfg["fc"],normalization=cfg["normalization"],
                 samples=records(qc),pca=records(pca),explained=explained,
                 complete_proteins=len(full),contrasts=cfg["comparisons"],
                 warnings=cfg.get("input_notes",[]),proteins=records(ann),statistics={},
                 abundance=json.loads(x.to_json(orient="index")),
                 input_scale=cfg["input_scale"])
    for co in cfg["comparisons"]:
        tab=pd.read_csv(r/"3.DEP"/(co["name"]+".all_proteins.tsv"),sep="\t")
        cols=["protein_id","log2FC","P.Value","adj.P.Val","direction","tested","n_control","n_treatment","n_complete_pairs"]
        payload["statistics"][co["name"]]=records(tab[cols])
    checks=enrich_payload(payload,r,cfg,resource_root)
    template=(ROOT/"report_template/report_interactive.html").read_text(encoding="utf-8")
    css=(ROOT/"report_template/report_interactive.css").read_text(encoding="utf-8")
    js=(ROOT/"report_template/report_interactive.js").read_text(encoding="utf-8")
    extension=(ROOT/"report_template/report_modules.js").read_text(encoding="utf-8")
    raw=json.dumps(payload,ensure_ascii=False,separators=(",",":"),allow_nan=False).replace("<","\\u003c")
    html=template.replace("/*INLINE_CSS*/",css).replace("/*INLINE_PLOTLY*/",Path(plotly_js).read_text(encoding="utf-8")).replace("/*INLINE_APP*/",js).replace("/*INLINE_MODULES*/",extension).replace("PAYLOAD_JSON",raw)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True);output.write_text(html,encoding="utf-8")
    manifest=dict(project=cfg["project"],protein_count=len(ann),sample_count=len(meta),
                  comparisons=[x["name"] for x in cfg["comparisons"]],
                  statistics="Existing limma results; presentation only",library="Plotly.js 3.1.0 (MIT)",
                  design_reference="https://alphapeptstats.readthedocs.io/en/latest/functions.html",
                  event_reference="https://plotly.com/javascript/click-events/",
                  source_report=str(report_dir),html_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                  plotly_sha256=hashlib.sha256(Path(plotly_js).read_bytes()).hexdigest(),
                  interaction_code="Custom linkage; AlphaPeptStats is not executed or bundled")
    manifest["validation"]=checks
    inputs=[r/"7.Provenance/config.resolved.json",r/"1.QC/samples.tsv",r/"1.QC/annotation.tsv",
            r/"5.Expression/normalized_log2.tsv",r/"5.Expression/sample_qc.tsv",r/"2.Cluster/PCA_scores.tsv"]
    inputs += [r/"3.DEP"/(co["name"]+".all_proteins.tsv") for co in cfg["comparisons"]]
    inputs += list((r/"6.Function").glob("*.tsv"))+list((r/"8.Pathways").glob("*.tsv"))
    manifest["input_sha256"]={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs if p.exists()}
    output.with_suffix(".manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    print(output)
if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--report-dir",required=True)
    ap.add_argument("--output",required=True);ap.add_argument("--plotly-js");ap.add_argument("--resource-root")
    a=ap.parse_args();build(a.report_dir,a.output,a.plotly_js,a.resource_root)
