#!/usr/bin/env python3
import sys,json,tarfile,re
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import proteomics as p
import pandas as pd,numpy as np
class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ["href","src"] and v and not v.startswith(("http:","https:","#","data:")):self.links.append(v)
reports=[]
for config in [p.ROOT/"examples/TMT_HepG2/config.v1.2.yaml",p.ROOT/"examples/DIA_cat/config.v1.2.yaml"]:
    c=p.load_config(config);o=Path(c["outdir"])
    v=json.loads((o/"0.InputQC/validation.json").read_text())
    stats=pd.read_csv(o/"2.Comparison/DEP_stat.tsv",sep="\t")
    x=pd.read_csv(o/"1.Preprocess/normalized_log2.tsv",sep="\t",index_col=0)
    m=pd.read_csv(o/"0.InputQC/samples.tsv",sep="\t",keep_default_na=False)
    state=json.loads((o/"log/state.json").read_text())
    assert set(state["stages"])==set(p.STAGES),state
    expected=(7608,6) if c["mode"]=="pd_tmt" else (607,9)
    assert (v["input_proteins"],v["samples"])==expected
    for co in c["comparisons"]:
        r=pd.read_csv(o/"2.Comparison"/(co["name"]+".all_proteins.tsv"),sep="\t")
        st=pd.read_csv(o/"2.Comparison"/(co["name"]+".statistics.tsv"),sep="\t")
        np.testing.assert_allclose(p.bh(st["P.Value"]),st["adj.P.Val"],rtol=1e-10,atol=1e-12)
        sub=m[m.group.isin([co["control"],co["treatment"]]) & (m.role=="sample")]
        if co["paired"]:
            deltas=[]
            for subject in sorted(sub.subject.unique()):
                aa=sub[(sub.subject==subject)&(sub.group==co["control"])]["sample"].iloc[0]
                bb=sub[(sub.subject==subject)&(sub.group==co["treatment"])]["sample"].iloc[0]
                deltas.append(x[bb]-x[aa])
            expected_fc=pd.concat(deltas,axis=1).mean(axis=1)
        else:
            expected_fc=x[sub[sub.group==co["treatment"]]["sample"]].mean(axis=1)-x[sub[sub.group==co["control"]]["sample"]].mean(axis=1)
        np.testing.assert_allclose(st.log2FC,expected_fc.loc[st.protein_id],atol=1e-9,rtol=1e-9)
        dep=pd.read_csv(o/"2.Comparison"/(co["name"]+".DEP.tsv"),sep="\t")
        assert (dep["adj.P.Val"]<=c["alpha"]).all()
        assert (dep.log2FC.abs()>=np.log2(c["fc"])).all()
        row=stats[stats.comparison==co["name"]].iloc[0]
        assert row.FDR_DEP==len(dep)==row.up+row.down
        assert r.protein_id.is_unique
    parser=Links();parser.feed((o/"proteomics_report.html").read_text())
    missing=[v for v in parser.links if not (o/unquote(v)).exists()]
    assert not missing,missing
    with tarfile.open(o/"proteomics_delivery.tar.gz") as t:
        assert "proteomics_report.html" in t.getnames()
        assert "proteomics_static_report.html" in t.getnames()
    doc=(o/"proteomics_report.html").read_text()
    payload=json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>',doc,re.S).group(1))
    for f in payload["extra"]["files"]+payload["extra"]["image_files"]:
        assert (o/f["path"]).exists(),f
    manifest=json.loads((o/"proteomics_report.manifest.json").read_text())
    assert manifest["validation"]["gsea_curves_verified"]>0
    old=Path(str(o).replace("pipeline_runs_v1.2_20261020","pipeline_runs_v1.1_20261020"))
    for rel in ["1.Preprocess/normalized_log2.tsv"]+[f"2.Comparison/{co['name']}.all_proteins.tsv" for co in c["comparisons"]]:
        assert p.sha(o/rel)==p.sha(old/rel),("Unexpected statistical change",rel)
    gx=pd.read_csv(o/"5.Pathways/GSVA_gene_matrix.tsv",sep="\t",index_col=0)
    assert np.isfinite(gx.to_numpy()).all() and gx.index.is_unique
    for co in c["comparisons"]:
        allp=pd.read_csv(o/"2.Comparison"/(co["name"]+".all_proteins.tsv"),sep="\t")
        ranked=allp[allp.tested & allp.gene.notna() & np.isfinite(allp.t)].groupby("gene").t.median()
        for gs in c["gene_sets"]:
            stem=co["name"]+"."+gs["name"]
            ranks=pd.read_csv(o/"5.Pathways"/(stem+".GSEA_ranking.tsv"),sep="\t",index_col=0)
            np.testing.assert_allclose(ranks.statistic,ranked.loc[ranks.index],atol=1e-10)
            gsea=o/"5.Pathways"/(stem+".GSEA.tsv")
            if gsea.exists():
                g=pd.read_csv(gsea,sep="\t")
                valid=g.pval.notna()
                np.testing.assert_allclose(p.bh(g.loc[valid,"pval"]),g.loc[valid,"padj"],atol=1e-10)
            scorefile=o/"5.Pathways"/(gs["name"]+".GSVA_scores.tsv")
            diff=o/"5.Pathways"/(stem+".GSVA_difference.tsv")
            if diff.exists():
                scores=pd.read_csv(scorefile,sep="\t",index_col=0)
                tab=pd.read_csv(diff,sep="\t",index_col=0)
                np.testing.assert_allclose(p.bh(tab["P.Value"]),tab["adj.P.Val"],atol=1e-10)
                sub=m[m.group.isin([co["control"],co["treatment"]]) & (m.role=="sample")]
                expected_diff=scores[sub[sub.group==co["treatment"]]["sample"]].mean(axis=1)-scores[sub[sub.group==co["control"]]["sample"]].mean(axis=1)
                np.testing.assert_allclose(tab.score_difference,expected_diff.loc[tab.index],atol=1e-10)
    reports.append({"project":c["project"],"input_proteins":v["input_proteins"],"samples":v["samples"],
                    "missing_report_links":missing,"contrasts":stats.to_dict(orient="records"),"status":"PASS"})
dest=p.ROOT/"validation/acceptance.v1.2.json"
p.jsave(dest,reports)
print(json.dumps(reports,ensure_ascii=False,indent=2))
