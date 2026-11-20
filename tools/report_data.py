"""Prepare existing QC, enrichment and pathway results for an offline report."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

def rows(path):
    p=Path(path)
    if not p.exists():return []
    return json.loads(pd.read_csv(p,sep="\t").to_json(orient="records",force_ascii=False))

def read_gmt(path):
    result={}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        a=line.split("\t")
        if len(a)>=3:result.setdefault(a[0],set()).update(a[2:])
    return result

def running_es(values,hits):
    n=len(values);nh=len(hits)
    if not 0<nh<n:raise ValueError("GSEA curve requires a nonempty proper subset")
    w=np.abs(np.asarray(values,float)[hits]);w=w/w.sum() if w.sum() else np.full(nh,1/nh)
    step=np.full(n,-1/(n-nh));step[hits]=w
    curve=np.cumsum(step)
    hi=max(0,float(curve.max()));lo=min(0,float(curve.min()))
    return hi if hi>-lo else lo

def enrich_payload(payload,r,cfg,resource_root=None):
    p=Path(r);extra=dict(ora={},gsea={},gsva={},ppi={},databases=[g["name"] for g in cfg["gene_sets"]])
    extra["summary"]=rows(p/"3.DEP/DEP_stat.tsv")
    extra["validation"]=json.loads((p/"1.QC/validation.json").read_text(encoding="utf-8"))
    extra["statuses"]={
        "function":json.loads((p/"6.Function/function_status.json").read_text()),
        "pathways":json.loads((p/"8.Pathways/pathway_status.json").read_text())}
    extra["correlation"]=rows(p/"2.Cluster/sample_correlation.tsv")
    qc=rows(p/"1.QC/protein_qc_metadata.tsv")
    extra["protein_qc"]=qc if qc and len(qc[0])>1 else []
    verified=0
    for gs in cfg["gene_sets"]:
        db=gs["name"];gmt=Path(gs["path"])
        if resource_root:gmt=Path(resource_root)/cfg["species"]/gmt.name
        sets=read_gmt(gmt) if gmt.exists() else {}
        extra["ora"][db]={};extra["gsea"][db]={}
        scorefile=p/"8.Pathways"/(db+".GSVA_scores.tsv")
        if scorefile.exists():
            table=pd.read_csv(scorefile,sep="\t",index_col=0)
            extra["gsva"][db]=dict(scores=json.loads(table.to_json(orient="index")),
                differences={co["name"]:rows(p/"8.Pathways"/(co["name"]+"."+db+".GSVA_difference.tsv")) for co in cfg["comparisons"]})
        for co in cfg["comparisons"]:
            name=co["name"];stem=name+"."+db
            extra["ora"][db][name]={direction:rows(p/"6.Function"/(stem+"."+direction+".tsv")) for direction in ["All","Up","Down"]}
            rr=p/"8.Pathways"/(stem+".GSEA_ranking.tsv")
            gr=p/"8.Pathways"/(stem+".GSEA.tsv")
            if not rr.exists() or not gr.exists():continue
            rank=pd.read_csv(rr,sep="\t");rankmap={g:i for i,g in enumerate(rank.gene)}
            vals=rank.statistic.to_numpy();result=rows(gr)
            for item in result:
                term=item["pathway"]
                if term not in sets:raise ValueError("Missing original gene set for GSEA curve: "+term)
                hits=sorted(rankmap[g] for g in sets[term] if g in rankmap)
                es=running_es(vals,hits)
                if item["ES"] is not None and not np.isclose(es,item["ES"],rtol=1e-6,atol=1e-7):
                    raise ValueError(f"GSEA ES mismatch: {name} {db} {term}: {es} != {item['ES']}")
                item["hits"]=hits;verified+=1
            extra["gsea"][db][name]=dict(rows=result,genes=rank.gene.tolist(),values=vals.tolist())
    for co in cfg["comparisons"]:
        edges=rows(p/"6.Function"/(co["name"]+".PPI.tsv"))
        graph=nx.Graph()
        for e in sorted(edges,key=lambda e:-e["score"])[:100]:graph.add_edge(e["gene1"],e["gene2"])
        layout=nx.spring_layout(graph,seed=42,weight=None,k=1.2) if len(graph)>2 else nx.circular_layout(graph)
        extra["ppi"][co["name"]]=dict(edges=edges,nodes=[dict(gene=g,x=float(layout[g][0]),y=float(layout[g][1]),degree=int(graph.degree(g))) for g in graph])
    extra["gsva_gene_count"]=len(rows(p/"8.Pathways/GSVA_gene_matrix.tsv"))
    extra["files"]=[dict(label=str(f.relative_to(p)),path="Report/"+f.relative_to(p).as_posix()) for f in sorted(p.rglob("*")) if f.is_file()]
    extra["image_files"]=[dict(label=f.name,path="image/"+f.name) for f in sorted(p.parent.joinpath("image").glob("*.pdf"))]
    payload["extra"]=extra
    return dict(gsea_curves_verified=verified,databases=extra["databases"])
