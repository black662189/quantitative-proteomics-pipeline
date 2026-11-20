"""Integration test: planted pathways, missing data, duplicate-gene collapse and direction."""
import json,sys,subprocess,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import proteomics as p
import numpy as np,pandas as pd
rng=np.random.default_rng(20261020)
with tempfile.TemporaryDirectory(prefix="proteomics_pathway_test_") as tmp:
    d=Path(tmp)
    for name in ["0.InputQC","1.Preprocess","2.Comparison"]:(d/name).mkdir()
    ids=["p"+str(i) for i in range(201)]
    genes=["G"+str(i) for i in range(200)]+["G0"]
    samples=["A1","A2","A3","B1","B2","B3"]
    x=rng.normal(10,1,(201,6));x[:20,3:]+=4;x[20:40,3:]-=4
    x[200]=x[0]+.5;x[199,0]=np.nan
    pd.DataFrame(x,index=pd.Index(ids,name="protein_id"),columns=samples).to_csv(d/"1.Preprocess/normalized_log2.tsv",sep="\t")
    pd.DataFrame(dict(protein_id=ids,gene=genes)).to_csv(d/"0.InputQC/annotation.tsv",sep="\t",index=False)
    pd.DataFrame(dict(sample=samples,group=["A"]*3+["B"]*3,subject=["1","2","3"]*2,
        role=["sample"]*6,plex=["p"]*6)).to_csv(d/"0.InputQC/samples.tsv",sep="\t",index=False)
    ranking=rng.normal(0,1,201);ranking[:20]+=8;ranking[20:40]-=8;ranking[200]=ranking[0]
    comparisons=[]
    for name,control,treat,sign in [("forward","A","B",1),("reverse","B","A",-1)]:
        comparisons.append(dict(name=name,control=control,treatment=treat,paired=True))
        pd.DataFrame(dict(protein_id=ids,gene=genes,t=ranking*sign,tested=True)).to_csv(
            d/"2.Comparison"/(name+".all_proteins.tsv"),sep="\t",index=False)
    gmt=d/"test.gmt"
    gmt.write_text("\n".join(name+"\tfixture\t"+"\t".join("G"+str(i) for i in indices)
        for name,indices in [("UP",range(20)),("DOWN",range(20,40)),("OTHER",range(40,80))]))
    cfg=dict(gsea=True,gsva=True,comparisons=comparisons,gene_sets=[dict(name="TEST",path=str(gmt),min_size=5,max_size=100)])
    (d/"config.json").write_text(json.dumps(cfg))
    subprocess.run(["/media/nfs/nfs02/yuzhij/miniconda3/envs/R453/bin/Rscript",str(p.ROOT/"scripts/pathways.R"),
                    str(d/"config.json"),str(d)],check=True)
    f=pd.read_csv(d/"5.Pathways/forward.TEST.GSEA.tsv",sep="\t",index_col=0)
    r=pd.read_csv(d/"5.Pathways/reverse.TEST.GSEA.tsv",sep="\t",index_col=0)
    assert f.loc["UP","NES"]>0 and f.loc["DOWN","NES"]<0
    assert r.loc["UP","NES"]<0 and r.loc["DOWN","NES"]>0
    assert f.loc["UP","padj"]<.05 and f.loc["DOWN","padj"]<.05
    a=pd.read_csv(d/"5.Pathways/forward.TEST.GSVA_difference.tsv",sep="\t",index_col=0)
    b=pd.read_csv(d/"5.Pathways/reverse.TEST.GSVA_difference.tsv",sep="\t",index_col=0)
    assert a.loc["UP","score_difference"]>0 and a.loc["DOWN","score_difference"]<0
    np.testing.assert_allclose(a.score_difference,-b.loc[a.index,"score_difference"],atol=1e-10)
    gx=pd.read_csv(d/"5.Pathways/GSVA_gene_matrix.tsv",sep="\t",index_col=0)
    assert "G199" not in gx.index and gx.index.is_unique
    np.testing.assert_allclose(gx.loc["G0"],np.median(x[[0,200]],axis=0))
    p.jsave(p.ROOT/"validation/pathway_synthetic.v1.1.json",
            dict(status="PASS",checks=["GSEA planted enrichment","GSEA reversed direction",
             "GSVA reversed paired effect","No missing-value imputation","Duplicate gene median"]))
    print("PASS: synthetic pathway integration checks")
