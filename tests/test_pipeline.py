import sys,unittest,tempfile,json,copy,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import proteomics as p
import numpy as np
import pandas as pd
class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.d=Path(self.tmp.name)
        self.df=pd.DataFrame({"protein_id":["p1","p2","p3"],"gene":["G1","G2;G3",""],
                             "a1":[1,2,0],"a2":[2,3,4],"b1":[4,5,6],"b2":[5,6,7]})
        self.m=pd.DataFrame({"sample":["A1","A2","B1","B2"],"column":["a1","a2","b1","b2"],
                             "group":["A","A","B","B"],"subject":["1","2","1","2"]})
        self.df.to_csv(self.d/"x.tsv",sep="\t",index=False);self.m.to_csv(self.d/"samples.csv",index=False)
        self.c={"_config":str(self.d/"config.yaml"),"input":{"path":"x.tsv"},"mode":"generic","samples":"samples.csv",
                "comparisons":[{"name":"A_vs_B","control":"A","treatment":"B","paired":True}],
                "input_scale":"linear","normalization":"none","min_replicates":2,"identification_q_max":.01,
                "min_unique_peptides":1,"input_notes":[],"upstream_normalization":"unknown"}
    def tearDown(self):self.tmp.cleanup()
    def test_bh(self):
        np.testing.assert_allclose(p.bh([.01,.04,.03,.8]),[.04,.053333333333,.053333333333,.8])
    def test_identifier_and_zero(self):
        p.input_stage(self.c,self.d);p.preprocess_stage(self.c,self.d)
        a=pd.read_csv(self.d/"0.InputQC/annotation.tsv",sep="\t",keep_default_na=False)
        self.assertEqual(a.gene.tolist(),["G1","",""])
        x=pd.read_csv(self.d/"1.Preprocess/normalized_log2.tsv",sep="\t",index_col=0)
        self.assertTrue(np.isnan(x.loc["p3","A1"]));self.assertEqual(x.loc["p1","B1"],2)
    def test_duplicate_protein_rejected(self):
        self.df.loc[1,"protein_id"]="p1";self.df.to_csv(self.d/"x.tsv",sep="\t",index=False)
        with self.assertRaisesRegex(ValueError,"duplicate"):p.input_stage(self.c,self.d)
    def test_unmatched_column_rejected(self):
        self.m.loc[0,"column"]="NOT_FOUND";self.m.to_csv(self.d/"samples.csv",index=False)
        with self.assertRaisesRegex(ValueError,"Unmatched"):p.input_stage(self.c,self.d)
    def test_unpaired_subject_rejected(self):
        self.m.loc[3,"subject"]="3";self.m.to_csv(self.d/"samples.csv",index=False)
        with self.assertRaisesRegex(ValueError,"pairing"):p.input_stage(self.c,self.d)
    def test_negative_abundance_rejected(self):
        self.df.loc[0,"a1"]=-1;self.df.to_csv(self.d/"x.tsv",sep="\t",index=False)
        with self.assertRaisesRegex(ValueError,"Negative"):p.input_stage(self.c,self.d)
    def test_taxonomy_guard(self):
        self.c.update({"taxid":9685,"gene_sets":[{"name":"GO","taxid":9606,"path":"missing.gmt"}]})
        with self.assertRaisesRegex(ValueError,"taxonomy"):p.function_stage(self.c,self.d)
    def test_limma_direction_and_pairing(self):
        rng=np.random.default_rng(12);n=120
        baseline=rng.normal(10,2,(n,1));subject=rng.normal(0,2,(n,3))
        a=baseline+subject+rng.normal(0,.05,(n,3));b=a+2+rng.normal(0,.05,(n,3))
        x=pd.DataFrame(np.c_[a,b],columns=["a1","a2","a3","b1","b2","b3"],index=["p"+str(i) for i in range(n)])
        x.index.name="protein_id";x.to_csv(self.d/"log.tsv",sep="\t")
        m=pd.DataFrame({"sample":x.columns,"group":["A"]*3+["B"]*3,"subject":["1","2","3"]*2,"role":["sample"]*6,"plex":["p"]*6})
        m.to_csv(self.d/"m.tsv",sep="\t",index=False)
        c={"comparisons":[{"name":"forward","control":"A","treatment":"B","paired":True},
                          {"name":"reverse","control":"B","treatment":"A","paired":True}],
           "min_replicates":2,"observed_fraction":.66}
        (self.d/"c.json").write_text(json.dumps(c))
        r="/media/nfs/nfs02/yuzhij/miniconda3/envs/metabolomics/bin/Rscript"
        if not Path(r).exists():self.skipTest("Server R required")
        subprocess.run([r,str(p.ROOT/"scripts/differential.R"),str(self.d/"c.json"),str(self.d/"log.tsv"),str(self.d/"m.tsv"),str(self.d/"stats")],check=True)
        f=pd.read_csv(self.d/"stats/forward.statistics.tsv",sep="\t");r=pd.read_csv(self.d/"stats/reverse.statistics.tsv",sep="\t")
        np.testing.assert_allclose(f.log2FC,-r.log2FC,atol=1e-9);self.assertTrue((f.log2FC>1.8).all())
        self.assertTrue((f["adj.P.Val"]<.05).all())
if __name__=="__main__":unittest.main(verbosity=2)
