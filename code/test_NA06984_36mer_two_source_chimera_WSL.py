#!/usr/bin/env python3
import gzip,csv
from pathlib import Path

BASE=Path("/mnt/c/IGOR/Papers/DELETIONS and ncRNAs")
OUT=BASE/"HGSV_183131_NA06984_chimera_forensics"; OUT.mkdir(parents=True,exist_ok=True)
HG=BASE/"HGSV_56081_extended_references"/"hg38.p14.fa.gz"
MT=BASE/"hg38_full_genome"/"chrM.fa.gz"
Q="CACTGTAAGCTAACTAGCATTAACCTTTTAAGTTAA"

def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def fasta(path):
    with gzip.open(path,"rt",errors="ignore") as f:
        name=None; seq=[]
        for line in f:
            if line.startswith(">"):
                if name is not None: yield name,"".join(seq).upper()
                name=line[1:].strip(); seq=[]
            else: seq.append(line.strip())
        if name is not None: yield name,"".join(seq).upper()
def locate(fragment, records):
    hits=[]
    for name,seq in records:
        for strand,x in [("+",fragment),("-",rc(fragment))]:
            st=0
            while True:
                i=seq.find(x,st)
                if i<0: break
                hits.append((name,i+1,i+len(x),strand))
                if len(hits)>=100: return hits
                st=i+1
    return hits

print("Loading references...",flush=True)
hg=list(fasta(HG)); mt=list(fasta(MT))
rows=[]
for cut in range(8,len(Q)-7):
    left,right=Q[:cut],Q[cut:]
    lh_hg=locate(left,hg); rh_hg=locate(right,hg)
    lh_mt=locate(left,mt); rh_mt=locate(right,mt)
    combos=[("nuclear","mt",lh_hg,rh_mt),("mt","nuclear",lh_mt,rh_hg),("nuclear","nuclear",lh_hg,rh_hg),("mt","mt",lh_mt,rh_mt)]
    for a,b,ha,hb in combos:
        if ha and hb:
            rows.append({"cut":cut,"left":left,"right":right,"left_source":a,"right_source":b,
                         "left_hits":len(ha),"right_hits":len(hb),"left_example":str(ha[0]),"right_example":str(hb[0])})
with open(OUT/"exact_two_source_decompositions.csv","w",newline="",encoding="utf-8") as f:
    fields=["cut","left","right","left_source","right_source","left_hits","right_hits","left_example","right_example"]
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
specific=[]
for r in rows:
    if "chr13" in r["left_example"] or "chr13" in r["right_example"]:
        if "mt" in (r["left_source"],r["right_source"]): specific.append(r)
with open(OUT/"summary.txt","w",encoding="utf-8") as f:
    f.write(f"Query: {Q}\n")
    f.write(f"Exact two-source decompositions (arms >=8 nt): {len(rows)}\n")
    f.write(f"chr13 + mt exact decompositions: {len(specific)}\n\n")
    for r in specific: f.write(str(r)+"\n")
print("Exact two-source decompositions:",len(rows))
print("chr13 + mt exact decompositions:",len(specific))
for r in specific[:20]: print(r)
print("DONE:",OUT)
