#!/usr/bin/env python3
import csv,gzip
from pathlib import Path
from collections import defaultdict,Counter

BASE=Path('/mnt/c/IGOR/Papers/DELETIONS and ncRNAs')
EXAMPLES=BASE/'HGSV_183131_all_GEUVADIS'/'HGSV_183131_examples.csv'
HG38=BASE/'HGSV_56081_extended_references'/'hg38.p14.fa.gz'
CHRM=BASE/'hg38_full_genome'/'chrM.fa.gz'
OUT=BASE/'HGSV_183131_chimera_across_GEUVADIS_fast'; OUT.mkdir(parents=True,exist_ok=True)
ANCHOR='AAGCTAACTAGCATTAACCT'
RNA36='CACTGTAAGCTAACTAGCATTAACCTTTTAAGTTAA'
KNOWN_LEFT='CACTGTAAGCTAACT'
KNOWN_RIGHT='AGCATTAACCTTTTAAGTTAA'
MINARM=8

def read_csv(p):
    with open(p,encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def rc(s): return s.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]

def fasta(path):
    with gzip.open(path,'rt',encoding='ascii',errors='ignore') as f:
        name=None; parts=[]
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name,''.join(parts).upper()
                name=line[1:].strip().split()[0]; parts=[]
            else: parts.append(line.strip())
        if name is not None: yield name,''.join(parts).upper()

rows=[]
for r in read_csv(EXAMPLES):
    ins=r.get('insert','').upper()
    if r.get('type') in {'exact36','anchor20_exact_insert','anchor20_contained'} or ANCHOR in ins or RNA36 in ins:
        rows.append(r)
print('Positive example rows to classify:',len(rows),flush=True)

frag_index=defaultdict(list); splits=defaultdict(list)
for rid,r in enumerate(rows):
    ins=r.get('insert','').upper()
    for cut in range(MINARM,len(ins)-MINARM+1):
        left,right=ins[:cut],ins[cut:]
        splits[rid].append((cut,left,right))
        for side,frag in (('L',left),('R',right)):
            frag_index[frag].append((rid,cut,side))
            frag_index[rc(frag)].append((rid,cut,side))

bylen=defaultdict(set)
for frag in frag_index: bylen[len(frag)].add(frag)

def scan(path,label):
    hits=defaultdict(list)
    for contig,seq in fasta(path):
        print(f'{label}: {contig}',flush=True)
        for L,wanted in bylen.items():
            if len(seq)<L: continue
            for i in range(len(seq)-L+1):
                mer=seq[i:i+L]
                if mer in wanted:
                    for rid,cut,side in frag_index[mer]:
                        k=(rid,cut,side)
                        if len(hits[k])<100:
                            hits[k].append((contig,i+1,i+L))
    return hits

nuc=scan(HG38,'nuclear')
mt=scan(CHRM,'mt')
classified=[]
for rid,r in enumerate(rows):
    ins=r.get('insert','').upper(); nucmt=[]; mtnuc=[]
    for cut,left,right in splits[rid]:
        ln=nuc.get((rid,cut,'L'),[]); rn=nuc.get((rid,cut,'R'),[])
        lm=mt.get((rid,cut,'L'),[]); rm=mt.get((rid,cut,'R'),[])
        if ln and rm: nucmt.append((cut,left,right,ln,rm))
        if lm and rn: mtnuc.append((cut,left,right,lm,rn))
    known=(ins==RNA36 or (ins.startswith(KNOWN_LEFT) and ins.endswith(KNOWN_RIGHT)))
    if known: verdict='EXACT_KNOWN_CHR17_MT_CHIMERA'
    elif nucmt: verdict='EXACT_NUCLEAR_TO_MT_CHIMERA'
    elif mtnuc: verdict='EXACT_MT_TO_NUCLEAR_CHIMERA'
    else: verdict='NO_EXACT_NUC_MT_SPLIT'
    best=max(nucmt,key=lambda x:min(len(x[1]),len(x[2]))) if nucmt else None
    classified.append({'sample':r.get('sample',''),'genotype_class':r.get('genotype_class',''),'insert':ins,
        'verdict':verdict,'n_nuclear_to_mt':len(nucmt),'n_mt_to_nuclear':len(mtnuc),
        'best_cut':best[0] if best else '','best_left':best[1] if best else '',
        'best_right':best[2] if best else '','best_left_example':str(best[3][0]) if best else '',
        'best_right_example':str(best[4][0]) if best else ''})
    print(f'[{rid+1}/{len(rows)}] {r.get("sample")} {r.get("genotype_class")} {verdict} nuc->mt={len(nucmt)}',flush=True)

with open(OUT/'HGSV_183131_positive_reads_chimera_classification.csv','w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(classified[0].keys())); w.writeheader(); w.writerows(classified)
stats=defaultdict(Counter)
for r in classified:
    g=r['genotype_class']; stats[g]['n']+=1
    if r['verdict']=='EXACT_KNOWN_CHR17_MT_CHIMERA': stats[g]['known']+=1
    if r['verdict'] in {'EXACT_KNOWN_CHR17_MT_CHIMERA','EXACT_NUCLEAR_TO_MT_CHIMERA'}: stats[g]['nucmt']+=1
    if r['verdict']=='NO_EXACT_NUC_MT_SPLIT': stats[g]['none']+=1
with open(OUT/'HGSV_183131_chimera_summary.txt','w',encoding='utf-8') as f:
    f.write(f'Positive example rows classified: {len(classified)}\n')
    f.write(f'Exact known chr17->chrM 36-mer: {sum(r["verdict"]=="EXACT_KNOWN_CHR17_MT_CHIMERA" for r in classified)}\n')
    f.write(f'Any exact nuclear->mt decomposition: {sum(r["verdict"] in {"EXACT_KNOWN_CHR17_MT_CHIMERA","EXACT_NUCLEAR_TO_MT_CHIMERA"} for r in classified)}\n')
    f.write(f'No exact nuclear->mt decomposition: {sum(r["verdict"]=="NO_EXACT_NUC_MT_SPLIT" for r in classified)}\n\n')
    for g,s in sorted(stats.items()):
        f.write(f'{g}: examples={s["n"]} known={s["known"]} nucmt={s["nucmt"]} none={s["none"]}\n')
print('DONE',OUT)
