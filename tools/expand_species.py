"""Reuse verified species-specific resources and fetch GO snapshots."""
import concurrent.futures, json, subprocess, sys, shutil, hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SPECS=[
 ("hsa",9606,"人","depends_hsa_GRCh38.p13_Ensembl_release_89"),
 ("fca",9685,"猫","depends_F.catus_Fca126_mat1.0_NCBI_release_105"),
 ("mmu",10090,"小鼠","depends_mmu_GRCm38.p6_Ensembl_release_96"),
 ("rno",10116,"大鼠","depends_rno_Rnor_6.0_Ensembl_release_90"),
 ("cfa",9615,"犬","depends_ROS_Cfam_1.0_NCBI_release_106"),
 ("ssc",9823,"猪","depends_Sscrofa11.1_NCBI_release_106"),
 ("mcf",9541,"食蟹猴","depends_Mfa_MFA1912RKSv2_NCBI_release_102"),
 ("mcc",9544,"恒河猴","depends_Mmul_10_NCBI_release_103")]
def build(spec):
 code,taxid,name,folder=spec;src=ROOT.parent/folder;dst=ROOT/"depends"/code;dst.mkdir(exist_ok=True)
 pairs=[(src/(code+".kegg_optimal.symbol.gmt"),dst/"KEGG.gmt")]
 ppi=list(src.glob("rename.uniq."+str(taxid)+".protein.links.*"))
 if code=="mcc":ppi=[src/"9544.protein_symbol.PPI.txt"]
 if len(ppi)!=1:raise ValueError(f"Ambiguous PPI resources: {code} {ppi}")
 pairs.append((ppi[0],dst/"PPI.tsv"))
 resources=[]
 for a,b in pairs:
  if not b.exists():shutil.copyfile(a,b)
  resources.append(dict(source=str(a),path=str(b),sha256=hashlib.file_digest(b.open("rb"),"sha256").hexdigest()))
 if not (dst/"GO_provenance.json").exists():
  subprocess.run([sys.executable,str(ROOT/"tools/download_go.py"),"--taxid",str(taxid),"--outdir",str(dst)],check=True)
 sets=[dict(name=n,taxid=taxid,path=str(dst/(n+".gmt")),min_size=5,max_size=500) for n in ["KEGG","GO_BP","GO_CC","GO_MF"]]
 import yaml
 config=dict(species=code,taxid=taxid,gene_sets=sets,ppi=dict(taxid=taxid,path=str(dst/"PPI.tsv")))
 (dst/"species.yaml").write_text(yaml.safe_dump(config,allow_unicode=True,sort_keys=False),encoding="utf-8")
 result=dict(code=code,taxid=taxid,name=name,resources=resources,config=str(dst/"species.yaml"),validation="database integrity only; biological example validated only for hsa/fca")
 (dst/"resource_provenance.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
 print("READY",code,flush=True);return result
if __name__=="__main__":
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(build,SPECS))
 (ROOT/"depends/species_registry.json").write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding="utf-8")
