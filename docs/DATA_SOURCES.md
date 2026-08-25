# Data sources

## GEUVADIS

RNA-seq data and sample metadata were drawn from the GEUVADIS lymphoblastoid-cell-line resource.

Key source publications:
- Lappalainen T, Sammeth M, Friedländer MR, et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature. 2013;501:506–511. doi:10.1038/nature12531.
- 't Hoen PAC, Friedländer MR, Almlöf J, et al. Reproducibility of high-throughput mRNA and small RNA sequencing across laboratories. Nature Biotechnology. 2013;31:1015–1022. doi:10.1038/nbt.2702.

## 1000 Genomes high-coverage WGS

High-coverage CRAM files were accessed from the 1000 Genomes/ENA public data collection and used for regional validation of NA06984 and NA19236.

## Human references expected by the scripts

The analysis used local copies of:
- hg38 / hg38.p14 nuclear reference sequence
- human mitochondrial sequence (chrM)
- GENCODE transcript reference
- EBV reference where applicable

Large public reference and raw sequencing files are intentionally **not included** in this repository.
