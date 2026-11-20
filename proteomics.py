#!/usr/bin/env python3
"""Reusable quantitative-matrix proteomics pipeline. Raw-MS search is not performed."""
import argparse, csv, hashlib, html, io, json, math, os, platform, re, shutil
import subprocess, sys, tarfile, time, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
import networkx as nx

ROOT=Path(__file__).resolve().parent
VERSION="1.2.0"
STAGES=["input","preprocess","qc","differential","function","pathways","report"]

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def jsave(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+".tmp");tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8");tmp.replace(p)

def tsv(p,df):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(p,sep="\t",index=False,na_rep="NA")

def bh(p):
    p=np.asarray(p,float);order=np.argsort(p);q=np.empty(len(p))
    if len(p):q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q

def pathvalue(c,v):
    p=Path(v)
    return p if p.is_absolute() else Path(c["_config"]).parent/p

def load_config(p):
    p=Path(p).resolve();c=yaml.safe_load(p.read_text(encoding="utf-8"));c["_config"]=str(p)
    if c.get("species_config"):
        species_path=pathvalue(c,c["species_config"])
        species=yaml.safe_load(species_path.read_text(encoding="utf-8"))
        for k in ["species","taxid","gene_sets","ppi"]:
            if k in species:c.setdefault(k,species[k])
    for k in ["project","species","taxid","mode","input","samples","comparisons","rscript","outdir"]:
        if k not in c:raise ValueError("Missing configuration: "+k)
    if c["mode"] not in ["pd_tmt","diann","generic"]:raise ValueError("mode must be pd_tmt, diann or generic")
    defaults={"input_scale":"linear","normalization":"none" if c["mode"]=="diann" else "median",
              "min_replicates":2,"observed_fraction":2/3,"fc":1.2,"alpha":0.05,
              "identification_q_max":0.01,"min_unique_peptides":1,"ppi_score":700,"gene_sets":[],
              "input_notes":[],"upstream_normalization":"unknown","gsea":True,"gsva":True,
              "pathway_rscript":"/media/nfs/nfs02/yuzhij/miniconda3/envs/R453/bin/Rscript"}
    for k,v in defaults.items():c.setdefault(k,v)
    if c["input_scale"] not in ["linear","log2"]:raise ValueError("input_scale must be linear or log2")
    if c["normalization"] not in ["none","median","bridge"]:raise ValueError("Invalid normalization")
    if not (c["fc"]>1 and 0<c["alpha"]<1 and 0<c["observed_fraction"]<=1 and c["min_replicates"]>=2):
        raise ValueError("Invalid statistical thresholds")
    names=[]
    for co in c["comparisons"]:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+",co["name"]):raise ValueError("Unsafe comparison name")
        co.setdefault("paired",False)
        if co["control"]==co["treatment"]:raise ValueError("Comparison groups must differ")
        names.append(co["name"])
    if len(names)!=len(set(names)) or not names:raise ValueError("Duplicate/empty comparisons")
    c["outdir"]=str(pathvalue(c,c["outdir"]).resolve())
    for gs in c["gene_sets"]:
        gs["path"]=str(pathvalue(c,gs["path"]).resolve())
        if str(gs["taxid"])!=str(c["taxid"]):raise ValueError("Gene set taxonomy mismatch")
    registry=ROOT/"depends/species_registry.json"
    if registry.exists():
        expected={r["code"]:r["taxid"] for r in json.loads(registry.read_text(encoding="utf-8"))}
        if c["species"] in expected and str(expected[c["species"]])!=str(c["taxid"]):
            raise ValueError("Species code / taxonomy mismatch")
    return c

def read_input(c):
    spec=c["input"];p=pathvalue(c,spec["path"])
    if p.suffix.lower()==".zip":
        with zipfile.ZipFile(p) as z:
            member=spec.get("member")
            if member not in z.namelist() or Path(member).name.startswith("~$"):raise ValueError("Explicit valid ZIP member required")
            info=z.getinfo(member)
            if info.file_size>1_000_000_000:raise ValueError("Matrix member >1 GB; extract and supply table explicitly")
            raw=io.BytesIO(z.read(member));ext=Path(member).suffix.lower()
    else:raw=p;ext=p.suffix.lower()
    if ext==".xlsx":
        if "sheet" not in spec:raise ValueError("XLSX requires explicit sheet; differential-only sheets are not valid input")
        return pd.read_excel(raw,sheet_name=spec["sheet"],engine="openpyxl")
    if ext not in [".tsv",".txt",".csv"]:raise ValueError("Supported input: .xlsx/.tsv/.csv (or selected ZIP member)")
    return pd.read_csv(raw,sep="," if ext==".csv" else "\t")

def input_stage(c,o):
    d=o/"0.InputQC";d.mkdir(parents=True,exist_ok=True)
    df=read_input(c);m=pd.read_csv(pathvalue(c,c["samples"]),dtype=str,keep_default_na=False)
    for k in ["sample","column","group"]:
        if k not in m:raise ValueError("samples.csv missing "+k)
    for k,v in {"subject":"","plex":"plex1","channel":"","role":"sample"}.items():
        if k not in m:m[k]=v
    if m["sample"].duplicated().any() or m["column"].duplicated().any():raise ValueError("Duplicate samples/columns")
    if not m["sample"].str.fullmatch(r"[A-Za-z0-9_.-]+").all():raise ValueError("Sample IDs must be filesystem-safe ASCII; retain original names in column")
    if not m.role.isin(["sample","reference","qc"]).all():raise ValueError("Unknown sample role")
    if any(not str(x).strip() for x in m["group"]):raise ValueError("Empty group")
    missing=set(m.column)-set(df.columns)
    if missing:raise ValueError("Unmatched sample columns: "+repr(missing))
    if c["mode"]=="pd_tmt":
        if not m.channel.str.len().gt(0).all() or m.duplicated(["plex","channel"]).any():raise ValueError("TMT channels must be unique within each plex")
        if m.plex.nunique()>1 and c["normalization"]!="bridge":raise ValueError("Multi-plex TMT requires bridge normalization and shared reference samples")
    if c["normalization"]=="bridge":
        if c["mode"]!="pd_tmt":raise ValueError("Bridge normalization is implemented for TMT only")
        for plex,g in m.groupby("plex"):
            if sum(g.role=="reference")!=1:raise ValueError("Require exactly one common reference per plex")
    bio=m[m.role=="sample"]
    for co in c["comparisons"]:
        sub=bio[bio.group.isin([co["control"],co["treatment"]])]
        if any(sum(sub.group==g)<c["min_replicates"] for g in [co["control"],co["treatment"]]):raise ValueError("Insufficient biological replicates: "+co["name"])
        if co["paired"]:
            if sub.plex.nunique()>1:raise ValueError("Paired multi-batch designs require a separately specified model; not supported in v1")
            if (sub.subject=="").any() or sub.duplicated(["subject","group"]).any():raise ValueError("Paired design requires one sample per subject per group")
            if set(sub[sub.group==co["control"]].subject)!=set(sub[sub.group==co["treatment"]].subject):raise ValueError("Incomplete subject pairing")
    aliases={"pd_tmt":("Accession","Gene Symbol","Description"),"diann":("Protein.Group","Genes","First.Protein.Description")}
    idcol,gcol,descol=aliases.get(c["mode"],(c["input"].get("id_column","protein_id"),c["input"].get("gene_column","gene"),c["input"].get("description_column","description")))
    for key,col in [("id_column",idcol),("gene_column",gcol),("description_column",descol)]:
        if key in c["input"]:
            if key=="id_column":idcol=c["input"][key]
            elif key=="gene_column":gcol=c["input"][key]
            else:descol=c["input"][key]
    if idcol not in df:raise ValueError("Protein ID column absent: "+idcol)
    ids=df[idcol].fillna("").astype(str)
    if (ids.str.strip()=="").any() or ids.duplicated().any():raise ValueError("Empty/duplicate protein IDs; resolve protein groups explicitly")
    numeric=df[list(m.column)].apply(pd.to_numeric,errors="coerce")
    invalid=df[list(m.column)].notna() & numeric.isna()
    if invalid.any().any() or np.isinf(numeric.to_numpy()).any():raise ValueError("Non-numeric/infinite abundance")
    if c["input_scale"]=="linear" and (numeric<0).any().any():raise ValueError("Negative linear abundance")
    genes=df[gcol].fillna("").astype(str) if gcol in df else pd.Series("",index=df.index)
    gene=genes.map(lambda s: s.strip() if s.strip() and not re.search(r"[;,| ]",s.strip()) else "")
    ann=pd.DataFrame({"protein_id":ids,"gene":gene,"gene_annotation":genes,
                      "description":df[descol].fillna("").astype(str) if descol in df else ""})
    # Protein groups remain independent units; ambiguous gene mappings are not expanded.
    keep=np.ones(len(df),dtype=bool);warnings=list(c["input_notes"]);filter_counts={}
    for col in ["Reverse","Potential contaminant","Only identified by site","Decoy","Contaminant"]:
        if col in df:
            bad=df[col].fillna("").astype(str).str.lower().isin(["+","true","1","yes"])
            filter_counts[col]=int(bad.sum());keep &= ~bad
    decoy=ids.str.match(r"^(REV__|CON__|DECOY_)",case=False)
    keep &= ~decoy;filter_counts["id_prefix"]=int(decoy.sum())
    qcol=c["input"].get("q_column","Exp. q-value: Combined" if c["mode"]=="pd_tmt" else "")
    if c["input"].get("q_column") and qcol not in df:raise ValueError("Configured identification q column absent")
    if qcol and qcol in df:
        q=pd.to_numeric(df[qcol],errors="coerce")
        if ((q.dropna()<0)|(q.dropna()>1)).any():raise ValueError("Invalid identification q value")
        ok=q.notna() & (q<=c["identification_q_max"]);filter_counts["identification_q"]=int((~ok).sum());keep &= ok
    else:warnings.append("输入矩阵没有可核验的鉴定q值；未重新控制质谱鉴定FDR，需核查上游搜库记录。")
    if "# Unique Peptides" in df:
        ok=pd.to_numeric(df["# Unique Peptides"],errors="coerce")>=c["min_unique_peptides"]
        filter_counts["unique_peptides"]=int((~ok).sum());keep &= ok
    if len(df)==0 or keep.sum()<2:raise ValueError("No usable protein matrix")
    numeric.columns=m["sample"].tolist()
    tsv(d/"abundance.tsv",pd.concat([ann[["protein_id"]],numeric],axis=1)[keep])
    tsv(d/"annotation.tsv",ann[keep]);tsv(d/"samples.tsv",m)
    # Quantitative protein-level QC only; peptide-level QC cannot be reconstructed from these inputs.
    qcmeta=[x for x in ["Coverage [%]","# Peptides","# PSMs","# Unique Peptides","MW [kDa]"] if x in df]
    tsv(d/"protein_qc_metadata.tsv",pd.concat([ann[["protein_id"]],df[qcmeta]],axis=1)[keep])
    warnings.extend(["当前输入为蛋白定量表；不包含原始谱图重搜库、肽段级质控或TMT同位素串扰校正。",
                     "样本分组及比较方向来自配置；正log2FC表示treatment相对control升高。",
                     "无标准样本/QC进样时，样本相关性不代表仪器稳定性验证。"])
    jsave(d/"validation.json",{"status":"PASS","input_proteins":len(df),"retained_proteins":int(keep.sum()),
          "samples":len(bio),"filter_counts":filter_counts,"ambiguous_or_missing_gene":int((ann.loc[keep,"gene"]=="").sum()),
          "warnings":warnings})

def preprocess_stage(c,o):
    d=o/"1.Preprocess";d.mkdir(exist_ok=True)
    x=pd.read_csv(o/"0.InputQC/abundance.tsv",sep="\t",index_col=0)
    m=pd.read_csv(o/"0.InputQC/samples.tsv",sep="\t",keep_default_na=False)
    if c["input_scale"]=="linear":x=np.log2(x.mask(x<=0))
    if not x.notna().any(axis=0).all():raise ValueError("At least one sample has no valid abundance")
    offsets={}
    if c["normalization"]=="median":
        # QC and bridge samples do not set the biological normalization target.
        for plex,mm in m.groupby("plex"):
            cols=mm[mm.role=="sample"]["sample"].tolist()
            if not cols:continue
            med=x[cols].median();target=med.median()
            for col in cols:offsets[col]=float(med[col]-target);x[col]-=offsets[col]
    elif c["normalization"]=="bridge":
        for plex,mm in m.groupby("plex"):
            ref=mm[mm.role=="reference"]["sample"].iloc[0]
            x[mm["sample"]]=x[mm["sample"]].subtract(x[ref].copy(),axis=0)
    x.index.name="protein_id";x.to_csv(d/"normalized_log2.tsv",sep="\t",na_rep="NA")
    count=pd.DataFrame({"sample":x.columns,"observed":x.notna().sum().values,"missing_fraction":x.isna().mean().values,
                        "median_log2":x.median().values})
    tsv(d/"sample_qc.tsv",m.merge(count,on="sample",validate="one_to_one"))
    jsave(d/"preprocessing.json",{"scale":"log2","normalization":c["normalization"],"offsets":offsets,
                                 "imputation_for_tests":"none","upstream_normalization":c["upstream_normalization"]})

def savefig(p):
    plt.tight_layout();plt.savefig(str(p)+".png",dpi=150,bbox_inches="tight");plt.savefig(str(p)+".pdf",bbox_inches="tight");plt.close()

def qc_stage(c,o):
    d=o/"3.Cluster";d.mkdir(exist_ok=True);image=o/"image";image.mkdir(exist_ok=True)
    x=pd.read_csv(o/"1.Preprocess/normalized_log2.tsv",sep="\t",index_col=0)
    m=pd.read_csv(o/"0.InputQC/samples.tsv",sep="\t",keep_default_na=False);m=m[m.role=="sample"]
    x=x[list(m["sample"])]
    plt.figure(figsize=(max(7,len(m)*.6),4));plt.boxplot([x[col].dropna() for col in x],tick_labels=list(x.columns),showfliers=False)
    plt.xticks(rotation=45,ha="right");plt.ylabel("log2 abundance");savefig(image/"abundance")
    plt.figure(figsize=(8,4));plt.bar(x.columns,x.isna().mean());plt.xticks(rotation=45,ha="right");plt.ylabel("Missing fraction");plt.ylim(0,1);savefig(image/"missing")
    corr=x.corr(min_periods=3);corr.to_csv(d/"sample_correlation.tsv",sep="\t",na_rep="NA")
    plt.figure(figsize=(7,6));plt.imshow(corr,vmin=-1,vmax=1,cmap="RdBu_r");plt.colorbar(label="Pearson r")
    plt.xticks(range(len(m)),m["sample"],rotation=90);plt.yticks(range(len(m)),m["sample"]);savefig(image/"correlation")
    full=x.dropna();full=full[full.var(axis=1)>0];status={"pca_complete_proteins":len(full)}
    if len(full)>=2:
        centered=full.to_numpy().T;centered-=centered.mean(axis=0)
        u,s,v=np.linalg.svd(centered,full_matrices=False);scores=u[:,:2]*s[:2];fraction=s*s/(s*s).sum()
        tab=pd.DataFrame(scores,columns=["PC1","PC2"]);tab.insert(0,"sample",m["sample"].values);tsv(d/"PCA_scores.tsv",tab)
        plt.figure(figsize=(7,5))
        for group in m.group.unique():
            take=(m.group==group).to_numpy();plt.scatter(scores[take,0],scores[take,1],label=group,s=60)
        for i,name in enumerate(m["sample"]):plt.annotate(name,(scores[i,0],scores[i,1]),fontsize=8)
        plt.xlabel(f"PC1 ({fraction[0]:.1%})");plt.ylabel(f"PC2 ({fraction[1]:.1%})");plt.legend();savefig(image/"pca")
        status["pca"]="complete cases only; centered, not variance-scaled"
    else:status["pca"]="SKIPPED: fewer than two complete variable proteins"
    jsave(d/"qc_status.json",status)
    metadata=pd.read_csv(o/"0.InputQC/protein_qc_metadata.tsv",sep="\t")
    for col,stem in [("Coverage [%]","coverage"),("# Unique Peptides","unique_peptides"),("MW [kDa]","molecular_weight")]:
        if col in metadata:
            vals=pd.to_numeric(metadata[col],errors="coerce").dropna()
            if len(vals):
                plt.figure(figsize=(7,4));plt.hist(vals,bins=40,color="#417AA0");plt.xlabel(col);plt.ylabel("Protein groups");savefig(image/stem)

def differential_stage(c,o):
    d=o/"2.Comparison";d.mkdir(exist_ok=True)
    rconfig=o/"log/resolved.json";jsave(rconfig,c)
    env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with open(o/"log/differential.log","w",encoding="utf-8") as log:
        subprocess.run([c["rscript"],str(ROOT/"scripts/differential.R"),str(rconfig),
                        str(o/"1.Preprocess/normalized_log2.tsv"),str(o/"0.InputQC/samples.tsv"),str(d)],
                       check=True,stdout=log,stderr=subprocess.STDOUT,env=env)
    annotation=pd.read_csv(o/"0.InputQC/annotation.tsv",sep="\t",keep_default_na=False)
    x=pd.read_csv(o/"1.Preprocess/normalized_log2.tsv",sep="\t",index_col=0)
    summary=[]
    for co in c["comparisons"]:
        name=co["name"];stats=pd.read_csv(d/(name+".statistics.tsv"),sep="\t")
        elig=pd.read_csv(d/(name+".eligibility.tsv"),sep="\t")
        r=annotation.merge(elig,on="protein_id",validate="one_to_one").merge(stats,on="protein_id",how="left",validate="one_to_one")
        r["FC"]=np.exp2(r.log2FC)
        effect=r.log2FC.abs()>=math.log2(c["fc"])
        r["direction"]=np.where(effect & (r["adj.P.Val"]<=c["alpha"]),np.where(r.log2FC>0,"Up","Down"),"NS")
        r.loc[~r.tested,"direction"]="Not_tested"
        r["exploratory_raw_p"]=effect & (r["P.Value"]<=c["alpha"])
        r=r.sort_values(["adj.P.Val","protein_id"],na_position="last")
        tsv(d/(name+".all_proteins.tsv"),r)
        sig=r[r.direction.isin(["Up","Down"])];tsv(d/(name+".DEP.tsv"),sig)
        tsv(d/(name+".exploratory_raw_p.tsv"),r[r.exploratory_raw_p])
        summary.append({"comparison":name,"control":co["control"],"treatment":co["treatment"],"paired":co["paired"],
                        "tested":int(r.tested.sum()),"up":int(sum(r.direction=="Up")),"down":int(sum(r.direction=="Down")),
                        "FDR_DEP":len(sig),"raw_p_candidates":int(r.exploratory_raw_p.sum())})
        plt.figure(figsize=(7,5))
        for label,color in [("NS","#A8B3C2"),("Up","#C74646"),("Down","#3269B0")]:
            a=r[r.direction==label];plt.scatter(a.log2FC,-np.log10(a["adj.P.Val"].clip(lower=1e-300)),s=10,c=color,label=label,alpha=.7)
        plt.axvline(math.log2(c["fc"]),ls="--",c="grey");plt.axvline(-math.log2(c["fc"]),ls="--",c="grey")
        plt.axhline(-math.log10(c["alpha"]),ls="--",c="grey");plt.xlabel("log2FC (treatment / control)");plt.ylabel("-log10 BH adjusted P");plt.legend();savefig(o/"image"/(name+".volcano"))
        if len(sig):
            ids=sig.protein_id.head(40);hm=x.loc[ids].copy();hm=hm.subtract(hm.mean(axis=1),axis=0).div(hm.std(axis=1).replace(0,np.nan),axis=0)
            plt.figure(figsize=(9,max(4,len(ids)*.22)));plt.imshow(np.ma.masked_invalid(hm),cmap="RdBu_r",vmin=-2,vmax=2,aspect="auto")
            plt.xticks(range(len(x.columns)),x.columns,rotation=90);plt.yticks(range(len(ids)),ids,fontsize=7);plt.colorbar(label="Row Z-score; grey = missing");plt.gca().set_facecolor("#dddddd");savefig(o/"image"/(name+".heatmap"))
    tsv(d/"DEP_stat.tsv",pd.DataFrame(summary))

def read_gmt(p):
    result={}
    for line in Path(p).read_text(encoding="utf-8").splitlines():
        a=line.split("\t")
        if len(a)>=3:
            if a[0] in result:result[a[0]]["genes"].update(a[2:])
            else:result[a[0]]={"description":a[1],"genes":set(a[2:])-{""}}
    return result

def function_stage(c,o):
    d=o/"4.Function";d.mkdir(exist_ok=True);statuses=[]
    for gs in c["gene_sets"]:
        if str(gs["taxid"])!=str(c["taxid"]):raise ValueError("Gene set taxonomy mismatch")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+",gs["name"]):raise ValueError("Invalid gene set name")
        sets=read_gmt(pathvalue(c,gs["path"]));mapped=set().union(*(x["genes"] for x in sets.values())) if sets else set()
        for co in c["comparisons"]:
            r=pd.read_csv(o/"2.Comparison"/(co["name"]+".all_proteins.tsv"),sep="\t",keep_default_na=False)
            tested=r[(r.tested==True)&(r.gene!="")]
            universe=set(tested.gene)&mapped
            terms={k:v for k,v in sets.items() if gs.get("min_size",5)<=len(v["genes"]&universe)<=gs.get("max_size",500)}
            for direction in ["All","Up","Down"]:
                take=tested.direction.isin(["Up","Down"]) if direction=="All" else tested.direction==direction
                hits=set(tested.loc[take,"gene"])&universe;rows=[]
                if hits and universe:
                    for term,v in terms.items():
                        bg=v["genes"]&universe;hit=bg&hits
                        rows.append({"term":term,"description":v["description"],"overlap":len(hit),"term_background":len(bg),
                                     "selected_genes":len(hits),"background_genes":len(universe),
                                     "pvalue":float(hypergeom.sf(len(hit)-1,len(universe),len(bg),len(hits))),"genes":";".join(sorted(hit))})
                res=pd.DataFrame(rows,columns=["term","description","overlap","term_background","selected_genes","background_genes","pvalue","genes"])
                res["padj"]=bh(res.pvalue);res=res.sort_values("padj")
                stem=co["name"]+"."+gs["name"]+"."+direction;tsv(d/(stem+".tsv"),res)
                statuses.append({"comparison":co["name"],"database":gs["name"],"direction":direction,"background":len(universe),
                                 "selected":len(hits),"tested_terms":len(res),"significant_terms":int(sum(res.padj<=c["alpha"])),
                                 "status":"DONE" if len(res) else "NO_TEST: no selected mapped genes or eligible terms"})
                show=res[res.overlap>0].head(15)
                if len(show):
                    plt.figure(figsize=(9,6));plt.barh(range(len(show)), -np.log10(show.padj.clip(lower=1e-300)),color="#417AA0")
                    plt.yticks(range(len(show)),[str(v)[:65] for v in show.term],fontsize=8);plt.gca().invert_yaxis()
                    plt.axvline(-math.log10(c["alpha"]),ls="--",c="grey");plt.xlabel("-log10 BH adjusted P (dashed = significance)");savefig(o/"image"/stem)
    if not c["gene_sets"]:statuses.append({"status":"SKIPPED: no gene-set database configured"})
    if c.get("ppi"):
        spec=c["ppi"]
        if str(spec["taxid"])!=str(c["taxid"]):raise ValueError("PPI taxonomy mismatch")
        needed=set()
        contrasts={}
        for co in c["comparisons"]:
            r=pd.read_csv(o/"2.Comparison"/(co["name"]+".DEP.tsv"),sep="\t",keep_default_na=False)
            contrasts[co["name"]]=set(r.gene)-{""};needed |= contrasts[co["name"]]
        edges={}
        with open(pathvalue(c,spec["path"]),encoding="utf-8") as f:
            for line in f:
                a=line.split()
                if len(a)!=3 or a[0] not in needed or a[1] not in needed or a[0]==a[1]:continue
                try:score=float(a[2])
                except ValueError:continue
                if score>=c["ppi_score"]:
                    pair=tuple(sorted(a[:2]));edges[pair]=max(score,edges.get(pair,0))
        for name,genes in contrasts.items():
            pairs=sorted([(a,b,s) for (a,b),s in edges.items() if a in genes and b in genes],key=lambda x:(-x[2],x[0],x[1]))
            tsv(d/(name+".PPI.tsv"),pd.DataFrame(pairs,columns=["gene1","gene2","score"]))
            if pairs:
                g=nx.Graph();g.add_weighted_edges_from(pairs[:100]);plt.figure(figsize=(8,7))
                pos=nx.circular_layout(g) if len(g)<=20 else nx.spring_layout(g,seed=42,weight=None,k=2/math.sqrt(len(g)))
                nx.draw_networkx(g,pos=pos,node_size=450,font_size=9,node_color="#97C6DD",edge_color="#6F8397",width=1.2)
                plt.axis("off");savefig(o/"image"/(name+".PPI"))
            statuses.append({"comparison":name,"database":"PPI","status":"DONE","edges":len(pairs)})
    else:statuses.append({"database":"PPI","status":"SKIPPED: no species-specific PPI database configured"})
    jsave(d/"function_status.json",statuses)

def pathways_stage(c,o):
    d=o/"5.Pathways";d.mkdir(exist_ok=True)
    if not c.get("gsea",True) and not c.get("gsva",True):
        jsave(d/"pathway_status.json",[{"status":"SKIPPED: disabled in config"}]);return
    resolved=o/"log/pathway_config.json";jsave(resolved,c)
    with open(o/"log/pathways.log","w",encoding="utf-8") as log:
        subprocess.run([c["pathway_rscript"],str(ROOT/"scripts/pathways.R"),str(resolved),str(o)],
                       check=True,stdout=log,stderr=subprocess.STDOUT,
                       env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1"))
    for p in d.glob("*.GSEA.tsv"):
        r=pd.read_csv(p,sep="\t").dropna(subset=["padj","NES"]).sort_values("padj").head(15)
        if len(r):
            plt.figure(figsize=(10,6));plt.barh(range(len(r)),r.NES,color=np.where(r.NES>0,"#C74646","#3269B0"))
            plt.yticks(range(len(r)),[str(v)[:65]+(" *" if q<=c["alpha"] else "") for v,q in zip(r.pathway,r.padj)],fontsize=8)
            plt.gca().invert_yaxis();plt.xlabel("NES; * BH adjusted P <= "+str(c["alpha"]));savefig(o/"image"/p.stem)
    for p in d.glob("*.GSVA_scores.tsv"):
        r=pd.read_csv(p,sep="\t",index_col=0)
        if len(r):
            r=r.loc[r.var(axis=1).nlargest(30).index]
            plt.figure(figsize=(10,8));plt.imshow(r,aspect="auto",cmap="RdBu_r")
            plt.xticks(range(len(r.columns)),r.columns,rotation=90);plt.yticks(range(len(r)),[str(v)[:65] for v in r.index],fontsize=7)
            plt.colorbar(label="GSVA score (top variance pathways)");savefig(o/"image"/p.stem)

def report_stage(c,o):
    report=o/"Report";report.mkdir(exist_ok=True)
    for src,dest in [("0.InputQC","1.QC"),("1.Preprocess","5.Expression"),("2.Comparison","3.DEP"),("3.Cluster","2.Cluster"),("4.Function","6.Function"),("5.Pathways","8.Pathways")]:
        shutil.copytree(o/src,report/dest,dirs_exist_ok=True)
    shutil.copytree(o/"log",report/"7.Provenance",dirs_exist_ok=True,ignore=shutil.ignore_patterns("run.lock"))
    shutil.copyfile(ROOT/"vendor/Plotly-LICENSE.txt",report/"7.Provenance/Plotly-LICENSE.txt")
    v=json.loads((o/"0.InputQC/validation.json").read_text());stat=pd.read_csv(o/"2.Comparison/DEP_stat.tsv",sep="\t")
    m=pd.read_csv(o/"0.InputQC/samples.tsv",sep="\t",keep_default_na=False)
    fn=json.loads((o/"4.Function/function_status.json").read_text())
    def esc(t):return html.escape(str(t))
    def table(df):
        return '<input class="search" placeholder="在此表中搜索" oninput="filterTable(this)">'+df.to_html(index=False,escape=True,border=0,classes="sortable")
    figure_counter=0
    def fig(stem,caption):
        nonlocal figure_counter
        path=o/"image"/(stem+".png")
        if not path.exists():return ""
        figure_counter+=1
        return '<figure><a href="image/'+esc(stem)+'.pdf"><img src="image/'+esc(stem)+'.png"></a><figcaption>图3-'+str(figure_counter)+' '+esc(caption)+'</figcaption></figure>'
    sections=[]
    sections.append(('summary','结果与结论',f'<p>本项目纳入 {v["samples"]} 个生物学样本，输入 {v["input_proteins"]} 个蛋白组，鉴定信息过滤后保留 {v["retained_proteins"]} 个。正式差异标准：BH校正P≤{c["alpha"]} 且FC≥{c["fc"]}或≤{1/c["fc"]:.4f}。各比较独立校正，不跨比较合并计数。</p>'+table(stat)+
                     '<p>raw_p_candidates 仅为未校正P值的探索性候选，不等同于FDR显著蛋白。零显著结果不会被替换为未校正结果。</p>'))
    sections.append(('overview','一、分析概述',f'<p>项目：{esc(c["project"])}；物种：{esc(c["species"])}（Taxonomy ID {esc(c["taxid"])}）；输入模式：{esc(c["mode"])}。</p><p>本流程从上游软件输出的蛋白丰度矩阵开始，用于相对定量和差异丰度分析。蛋白组保留原始ID，含多个或缺失基因符号的条目不强行展开用于基因集富集。</p>'+table(m[["sample","group","subject","plex","channel","role"]])))
    sections.append(('workflow','二、工作流程','<p>输入与样本核验 → 鉴定信息过滤 → log2转换及配置指定归一化 → 缺失率、相关性和PCA质控 → 各比较独立缺失过滤及limma建模 → BH校正 → 功能富集与PPI → 结果整理与报告。</p><p>模式适配只负责读取定量矩阵，不再次执行TMT报告离子提取、DIA谱图搜索或肽段到蛋白汇总。</p>'))
    results='<h3>3.1 数据质控</h3>'+table(pd.read_csv(o/"1.Preprocess/sample_qc.tsv",sep="\t"))+fig('abundance','归一化后log2丰度分布；箱线图仅使用观测值。')+fig('missing','每个样本的缺失比例；线性输入中的零值视为未检出。')+fig('correlation','生物学样本Pearson相关性；使用成对共同观测蛋白，少于3个共同观测时为空。')+fig('pca','PCA仅使用所有生物学样本均有观测且方差非零的蛋白；未用填补值制造分离。')
    results+=fig("coverage","蛋白序列覆盖率分布，来自上游鉴定表。")+fig("unique_peptides","每个蛋白的唯一肽段数分布，来自上游鉴定表。")+fig("molecular_weight","蛋白分子量分布，来自上游鉴定表。")
    results+='<p>PCA处理记录：'+esc((o/"3.Cluster/qc_status.json").read_text())+'</p>'
    for co in c["comparisons"]:
        name=co["name"];r=pd.read_csv(o/"2.Comparison"/(name+".DEP.tsv"),sep="\t")
        results+=f'<h3 id="{esc(name)}">3.2 {esc(name)}</h3><p>方向：{esc(co["treatment"])}/{esc(co["control"])}；设计：'+('配对（处理减对照的配对差值模型）' if co["paired"] else '独立组（多批次时纳入plex协变量）')+'</p>'
        results+=fig(name+'.volcano','火山图：横轴为log2FC，纵轴为−log10 BH校正P。红色升高，蓝色降低，灰色未达正式差异标准。')
        results+=f'<p>正式差异蛋白 {len(r)} 个。<a href="Report/3.DEP/{esc(name)}.all_proteins.tsv">下载全部蛋白结果</a>；<a href="Report/3.DEP/{esc(name)}.DEP.tsv">下载正式差异结果</a>。</p>'
        if len(r):results+=table(r[["protein_id","gene","log2FC","P.Value","adj.P.Val","direction"]].head(100))+fig(name+'.heatmap','最多展示40个正式差异蛋白的行Z-score，灰色为缺失，不用于检验。')
        for gs in c["gene_sets"]:
            results+=fig(name+'.'+gs["name"]+'.All',gs["name"]+' 富集前15项，虚线为BH阈值；展示条目不一定显著。')
        results+=fig(name+'.PPI','正式差异基因间的高置信相互作用，图中最多展示100条边；全量边表另存。')
    results+='<h3>3.3 功能分析执行记录</h3>'+table(pd.DataFrame(fn).fillna(""))
    results+='<h3>3.4 GSEA与GSVA</h3><p>GSEA使用全部可检验蛋白的有方向limma t统计量，同一基因取中位数；不按差异阈值预筛。正NES表示处理组方向富集。GSVA将同一基因的蛋白log2丰度按样本取中位数，仅保留所有生物学样本均有观测的非恒定基因，不填补；高缺失数据的通路覆盖因此可能减少。GSVA分数差异是处理减对照，不是蛋白Fold Change。</p>'
    results+=table(pd.DataFrame(json.loads((o/"5.Pathways/pathway_status.json").read_text())).fillna(""))
    for p in sorted((o/"5.Pathways").glob("*.GSEA.tsv")):
        r=pd.read_csv(p,sep="\t")
        results+='<h4>'+esc(p.stem)+'</h4>'+table(r.head(15))+fig(p.stem,'按BH校正P排序的前15条通路；星号代表达到设定显著性阈值，其余仅展示。')
    for p in sorted((o/"5.Pathways").glob("*.GSVA_scores.tsv")):
        results+=fig(p.stem,'GSVA分数热图：展示跨样本方差最大的30条通路，不等同于显著差异通路。')
    for p in sorted((o/"5.Pathways").glob("*.GSVA_difference.tsv")):
        results+='<h4>'+esc(p.stem)+'</h4>'+table(pd.read_csv(p,sep="\t").head(15))
    sections.append(('results','三、分析结果',results))
    methods=f'<p>输入尺度：{esc(c["input_scale"])}；归一化：{esc(c["normalization"])}；上游归一化记录：{esc(c["upstream_normalization"])}。线性丰度先将零视为缺失后log2转换，统计检验不填补缺失值。</p><p>每个比较中，两组均须至少{c["min_replicates"]}个且至少{c["observed_fraction"]:.0%}样本有观测。配对比较另要求相同比例的完整配对。配对检验对逐个体处理减对照差值拟合截距，独立比较拟合组别及可识别的批次效应。limma经验贝叶斯方差收缩使用trend=TRUE、robust=FALSE；每个比较在全部受检蛋白内执行BH校正。</p><p>富集使用超几何检验；背景为本比较受检蛋白中具有唯一基因符号且可映射到当前数据库的基因集合。对符合背景大小条件的全部条目（含零命中项）进行BH校正，分别计算全部、上调、下调集合。未注释蛋白及歧义映射不进入富集。</p>'
    methods+='<p>示例GO库使用UniProt直接注释，不传播祖先节点；人源库限reviewed蛋白，猫库包含全部条目。KEGG为现有转录组物种GMT快照，原始下载日期未知；PPI来源和版本见流程数据库来源记录。</p><p>方法参考：<a href="https://bioconductor.org/packages/limma/">limma</a>；<a href="https://docs.quantms.org/workflows/">quantms定量流程分支</a>。本流程实现矩阵下游分析，不声称运行了quantms原始谱图分析。</p>'
    methods+='<p>GSEA实际执行fgseaMultilevel，固定随机种子20261020，eps=0，以limma t统计量排序，BH按每个比较/数据库单独校正。GSVA实际执行GSVA::gsvaParam及gsva，Gaussian核；通路分数用limma建模，配对和批次设计与蛋白分析一致，trend=FALSE，按每个比较/数据库BH校正。通路大小按实际可映射基因过滤。完整基因矩阵、排名、leading-edge及软件版本随报告归档。</p>'
    sections.append(('methods','四、方法描述',methods))
    sections.append(('files','五、结果文件说明','<p>表格为UTF-8 TSV，保留真实扩展名，可用Excel导入。点击表头可排序，搜索框仅筛选当前显示表格；报告差异表最多内嵌100行，全量请下载。</p><ul>'+''.join('<li><a href="'+esc(str(p.relative_to(o)).replace(os.sep,'/'))+'">'+esc(str(p.relative_to(report)))+'</a></li>' for p in sorted(report.rglob('*')) if p.is_file())+'</ul>'))
    sections.append(('limits','六、局限与建议','<ul>'+''.join('<li>'+esc(t)+'</li>' for t in v["warnings"])+
                     '</ul><p>建议对候选蛋白做独立验证。报告仅陈述实际计算结果，不将PCA分离、未校正P值候选或通路展示等同于已证实生物学机制。</p>'))
    css=(ROOT/"report_template/style.css").read_text(encoding="utf-8")
    script=(ROOT/"report_template/table.js").read_text(encoding="utf-8")
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+esc(c["project"])+'</title><style>'+css+'</style><aside><h2>分析报告目录</h2>'+''.join('<a href="#'+sid+'">'+title+'</a>' for sid,title,_ in sections)+'</aside><main><header><h1>'+esc(c["project"])+'</h1><p>定量蛋白质组分析报告 · 流程版本 '+VERSION+'</p></header>'+''.join('<section id="'+sid+'"><h2>'+title+'</h2>'+body+'</section>' for sid,title,body in sections)+'</main><script>'+script+'</script></html>'
    (o/"proteomics_static_report.html").write_text(doc,encoding="utf-8")
    from tools.build_interactive_report import build
    build(o,o/"proteomics_report.html")
    with tarfile.open(o/"proteomics_delivery.tar.gz","w:gz") as tar:
        for n in ["proteomics_report.html","proteomics_static_report.html","proteomics_report.manifest.json","Report","image"]:tar.add(o/n,arcname=n)

def provenance(c):
    files=[Path(c["_config"]),pathvalue(c,c["samples"]),pathvalue(c,c["input"]["path"]),ROOT/"proteomics.py",ROOT/"scripts/differential.R",ROOT/"scripts/pathways.R",ROOT/"report_template/style.css",ROOT/"report_template/table.js"]
    files += [ROOT/p for p in ["tools/build_interactive_report.py","tools/report_data.py","report_template/report_interactive.html",
                "report_template/report_interactive.css","report_template/report_interactive.js","report_template/report_modules.js",
                "vendor/plotly-3.1.0.min.js","vendor/Plotly-LICENSE.txt"]]
    if c.get("species_config"):files.append(pathvalue(c,c["species_config"]))
    files += [pathvalue(c,x["path"]) for x in c["gene_sets"]]
    if c.get("ppi"):files.append(pathvalue(c,c["ppi"]["path"]))
    hashes={str(p.resolve()):sha(p) for p in files}
    return {"version":VERSION,"python":sys.version,"platform":platform.platform(),"input_and_code_sha256":hashes,
            "packages":{p:__import__(p).__version__ for p in ["numpy","pandas","scipy","matplotlib","yaml","networkx","openpyxl"]}}

def run(c,force=False,only=None):
    o=Path(c["outdir"])
    if o.exists() and any(o.iterdir()) and not (o/"log/state.json").exists():
        raise ValueError("Nonempty output directory without pipeline state; choose a new output directory")
    o.mkdir(parents=True,exist_ok=True);(o/"log").mkdir(exist_ok=True)
    import fcntl
    with open(o/"log/run.lock","w") as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError("Another pipeline process owns this output directory")
        prov=provenance(c);fingerprint=hashlib.sha256(json.dumps([c,prov],sort_keys=True).encode()).hexdigest()
        statepath=o/"log/state.json";state=json.loads(statepath.read_text()) if statepath.exists() else {}
        if state.get("fingerprint")!=fingerprint:state={"fingerprint":fingerprint,"stages":{}}
        jsave(statepath,state)
        jsave(o/"log/provenance.json",prov);jsave(o/"log/config.resolved.json",c)
        if only=="report" and not all(x in state.get("stages",{}) for x in STAGES[:-1]):raise ValueError("Report refresh requires valid prior analysis state")
        for stage in STAGES:
            if only and stage!=only:continue
            prev=state["stages"].get(stage)
            if not force and not only and prev and all(Path(p).exists() and sha(p)==h for p,h in prev["outputs"].items()):
                print("RESUME",stage,flush=True);continue
            # Invalidated stage invalidates successors, including old report success.
            for later in STAGES[STAGES.index(stage):]:state["stages"].pop(later,None)
            jsave(statepath,state)
            # Never leave a stale success report available after an invalidated/failed analysis.
            for old in [o/"proteomics_report.html",o/"proteomics_static_report.html",o/"proteomics_report.manifest.json",o/"proteomics_delivery.tar.gz"]:
                if old.exists():old.unlink()
            owned={"input":"0.InputQC","preprocess":"1.Preprocess","qc":"3.Cluster",
                   "differential":"2.Comparison","function":"4.Function","pathways":"5.Pathways","report":"Report"}
            target=o/owned[stage]
            if target.exists():shutil.rmtree(target)
            qc_stems={"abundance","missing","correlation","pca","coverage","unique_peptides","molecular_weight"}
            def image_owner(f):
                if f.stem.endswith((".GSEA",".GSVA_scores")):return "pathways"
                return "qc" if f.stem in qc_stems else ("differential" if f.stem.endswith((".volcano",".heatmap")) else "function")
            if stage in ["qc","differential","function","pathways"] and (o/"image").exists():
                for f in (o/"image").iterdir():
                    owner=image_owner(f)
                    if f.is_file() and owner==stage:f.unlink()
            print("RUN",stage,flush=True);started=time.time()
            globals()[stage+"_stage"](c,o)
            targets={"input":["0.InputQC"],"preprocess":["1.Preprocess"],"qc":["3.Cluster","image"],
                     "differential":["2.Comparison","image"],"function":["4.Function","image"],
                     "pathways":["5.Pathways","image"],"report":["Report","proteomics_report.html","proteomics_static_report.html","proteomics_report.manifest.json","proteomics_delivery.tar.gz"]}[stage]
            outputs={}
            for target in targets:
                p=o/target
                for f in ([p] if p.is_file() else p.rglob("*")):
                    # Images tracked by producing stages only, avoiding downstream cache invalidation.
                    if f.is_file() and (target!="image" or image_owner(f)==stage):
                        outputs[str(f)]=sha(f)
            state["stages"][stage]={"seconds":round(time.time()-started,3),"outputs":outputs};jsave(statepath,state)
        print("COMPLETE",o/"proteomics_report.html",flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument("command",choices=["run","validate","report"]);a.add_argument("-c","--config",required=True);a.add_argument("--force",action="store_true")
    args=a.parse_args();c=load_config(args.config)
    if args.command=="validate":
        import tempfile
        with tempfile.TemporaryDirectory(prefix="proteomics_validate_") as tmp:
            input_stage(c,Path(tmp));print((Path(tmp)/"0.InputQC/validation.json").read_text())
    else:run(c,args.force,"report" if args.command=="report" else None)
if __name__=="__main__":
    try:main()
    except Exception as e:
        print("FAILED:",e,file=sys.stderr);raise
