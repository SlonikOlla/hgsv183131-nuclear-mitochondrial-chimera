# Analysis workflow

1. **Genotype-stratified RNA screen**
   - Screen GEUVADIS libraries for the HGSV_183131 36-mer and diagnostic 20-mer.
   - Classify donors as stringent zero-copy Tier A, other zero-copy, one-copy, two-copy, or missing/other.
   - Script: `screen_HGSV_183131_all_GEUVADIS_WSL.py`.

2. **WGS validation**
   - Validate RNA-positive zero-copy donors using high-coverage 1000 Genomes WGS.
   - Measure deletion depth, flanking depth, zero-depth fraction, and exact 20-mer coverage.
   - Script: `validate_HGSV_183131_two_RNA_positive_donors_WGS_v2.py`.

3. **Single-source sequence-origin forensics**
   - Search recurrent 36-mer and 20-mer against extended human sequence, transcriptome, EBV, and technical motifs.
   - Script: `forensic_HGSV_183131_sequence_origin_WSL.py`.

4. **Human mitochondrial + raw-read audit**
   - Compare the recurrent 36-mer base-by-base against human chrM.
   - Inspect raw-read quality at RNA-vs-mtDNA mismatch positions.
   - Script: `forensic_NA06984_human_mtDNA_and_raw_reads_WSL.py`.

5. **Exact two-source decomposition**
   - Test every split with both arms >=8 nt for nuclear→mt, mt→nuclear, nuclear→nuclear, and mt→mt explanations.
   - Script: `test_NA06984_36mer_two_source_chimera_WSL.py`.

6. **Population-scale chimera audit**
   - Apply exact two-source reconstruction to all retained HGSV_183131-positive examples.
   - Final memory-efficient script: `audit_HGSV_183131_positive_reads_for_chimera_WSL_fast.py`.

## Stopping rule

The RNA-cache interpretation was rejected for HGSV_183131 once:
- WGS confirmed the genomic deletion,
- the apparent RNA signal remained detectable,
- and all positive stringent zero-copy examples were exactly explainable by nuclear–mitochondrial sequence architecture.

The analysis therefore distinguishes **valid RNA-seq observations** from **incorrect genomic provenance assignment**.
