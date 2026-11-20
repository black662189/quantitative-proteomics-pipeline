#!/usr/bin/env python3
"""Build versioned species-specific direct GO annotation GMT snapshots from UniProt."""
import urllib.request,urllib.parse,csv,io,re,json,hashlib,datetime,argparse
from pathlib import Path
a=argparse.ArgumentParser();a.add_argument('--taxid',required=True);a.add_argument('--outdir',required=True)
args=a.parse_args();d=Path(args.outdir);d.mkdir(parents=True,exist_ok=True)
query='organism_id:'+args.taxid
if args.taxid=='9606':query+=' AND reviewed:true'
url='https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode({'query':query,'format':'tsv','fields':'accession,gene_primary,go_p,go_c,go_f'})
with urllib.request.urlopen(url,timeout=120) as r:raw=r.read();headers=dict(r.headers)
(d/'uniprot_GO_source.tsv').write_bytes(raw)
rows=list(csv.DictReader(io.StringIO(raw.decode()),delimiter='\t'));counts={}
for kind,col in [('GO_BP','Gene Ontology (biological process)'),('GO_CC','Gene Ontology (cellular component)'),('GO_MF','Gene Ontology (molecular function)')]:
 terms={}
 for row in rows:
  gene=row['Gene Names (primary)'].strip()
  if not gene or re.search(r'[ ;,|]',gene):continue
  for label,go in re.findall(r'([^;]+?) \[(GO:\d+)\]',row.get(col,'')):
   term=go+' '+label.strip();terms.setdefault(term,set()).add(gene)
 with (d/(kind+'.gmt')).open('w',encoding='utf-8') as f:
  for term,genes in sorted(terms.items()):f.write(term+'\tUniProt_direct_annotation\t'+'\t'.join(sorted(genes))+'\n')
 counts[kind]=len(terms)
(d/'GO_provenance.json').write_text(json.dumps({'url':url,'taxid':args.taxid,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(raw).hexdigest(),'headers':headers,'counts':counts,'scope':'Direct UniProt annotations, no ancestor propagation; human reviewed proteins only, other species all entries.'},ensure_ascii=False,indent=2))
print(counts,flush=True)
