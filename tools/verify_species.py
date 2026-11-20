import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import proteomics as p
import yaml
records=json.loads((p.ROOT/"depends/species_registry.json").read_text())
assert len(records)==8 and len({r["taxid"] for r in records})==8
checks=[]
for r in records:
    c=yaml.safe_load(Path(r["config"]).read_text())
    assert c["species"]==r["code"] and c["taxid"]==r["taxid"]
    counts={};genes=set()
    for gs in c["gene_sets"]:
        assert gs["taxid"]==r["taxid"]
        sets=p.read_gmt(gs["path"])
        assert len(sets)>0
        counts[gs["name"]]=len(sets)
        if gs["name"]=="KEGG":genes=set().union(*(v["genes"] for v in sets.values()))
    go=json.loads((Path(r["config"]).parent/"GO_provenance.json").read_text())
    assert str(go["taxid"])==str(r["taxid"])
    seen=set();edges=0
    with open(c["ppi"]["path"]) as f:
        for line in f:
            a=line.split()
            if len(a)!=3:raise ValueError((r["code"],"Invalid PPI",line[:100]))
            try:score=float(a[2])
            except ValueError:
                if edges==0:continue
                raise
            assert 0<=score<=1000
            seen.update(a[:2]);edges+=1
    assert edges>0 and len(seen&genes)>0
    checks.append(dict(species=r["code"],taxid=r["taxid"],terms=counts,
        ppi_edges=edges,kegg_ppi_shared_symbols=len(seen&genes),status="PASS"))
p.jsave(p.ROOT/"validation/species.v1.1.json",checks)
print(json.dumps(checks,ensure_ascii=False,indent=2))
