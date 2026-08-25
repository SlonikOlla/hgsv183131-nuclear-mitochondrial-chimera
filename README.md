# Recurrent nuclear–mitochondrial chimeric RNA-seq reads can mimic transcription from a homozygously deleted human genomic interval

This repository contains the analysis code, derived data, figures, and manuscript-supporting files for a re-analysis of **HGSV_183131** in GEUVADIS/1000 Genomes data.

## Scientific background

A previous study reported non-coding RNA reads matching homozygously deleted genomic regions and proposed that residual RNA sequence information might persist after loss of the corresponding DNA. The present analysis revisits that interpretation with stronger WGS validation and sequence-origin forensics.

For HGSV_183131, the deletion is independently supported by WGS in RNA-positive zero-copy donors, including zero genomic coverage across the diagnostic 20-mer. However, the recurrent RNA sequence can be reconstructed exactly as a **nuclear–mitochondrial chimera**. Across 83 retained HGSV_183131-positive RNA examples, 77 have an exact nuclear-to-mitochondrial decomposition, and all five positive examples in the stringent zero-copy Tier-A group have such an explanation.

The data therefore do **not** support an RNA-cache interpretation for HGSV_183131. They demonstrate a sequence-provenance confounder: a chimeric junction can recreate a short k-mer that exactly matches an unrelated genomic locus.

## Repository layout

- `code/` — analysis scripts used for genotype-stratified RNA screening, WGS validation, sequence-origin forensics, mitochondrial comparison, exact two-source decomposition, and the final 83-read chimera audit.
- `data/derived/` — final classification tables, QC table, and WGS validation summary.
- `data/forensics/` — derived sequence-origin, mitochondrial, and two-source forensic outputs.
- `figures/` — manuscript figures and an editable PowerPoint copy of Figure 2.
- `paper/` — current manuscript draft, Supplementary Methods, and Supplementary Table S1.
- `docs/ANALYSIS_WORKFLOW.md` — step-by-step analysis logic and interpretation.
- `docs/DATA_SOURCES.md` — public data sources and local reference expectations.

## Key sequences

- Recurrent 36-mer: `CACTGTAAGCTAACTAGCATTAACCTTTTAAGTTAA`
- HGSV_183131 diagnostic 20-mer: `AAGCTAACTAGCATTAACCT`
- Exact recurrent decomposition:
  - nuclear arm: `CACTGTAAGCTAACT`
  - mitochondrial arm: `AGCATTAACCTTTTAAGTTAA`

## Key final result

See `data/derived/HGSV_183131_positive_reads_chimera_classification.csv`.

Final audit:
- positive examples classified: 83
- exact nuclear-to-mitochondrial reconstruction: 77/83 (92.8%)
- recurrent known chr17-to-chrM 36-mer: 56/83 (67.5%)
- positive stringent zero-copy Tier-A examples with exact nuclear-to-mitochondrial reconstruction: 5/5

## Reproducibility notes

The scripts were developed and run under WSL/Linux and assume `samtools`, Python 3, and access to public ENA/1000 Genomes resources. Several scripts contain project-specific absolute paths under:

`/mnt/c/IGOR/Papers/DELETIONS and ncRNAs`

Before reuse on another system, edit the `BASE` variable near the top of each script or preserve the same directory layout.

The final memory-efficient chimera audit is:

`code/audit_HGSV_183131_positive_reads_for_chimera_WSL_fast.py`

## Public data

Primary source datasets are GEUVADIS and 1000 Genomes resources. See `docs/DATA_SOURCES.md` and the manuscript references for source publications.

## Interpretation

Exact two-source decomposition establishes a plausible sequence architecture and invalidates assignment of the junction-spanning 20-mer specifically to HGSV_183131. It does **not** determine whether the recurrent nuclear–mitochondrial chimeras were present biologically in the cells or generated during reverse transcription/library preparation.

## Citation

Please cite the associated manuscript once published. A `CITATION.cff` file is included for repository citation before publication.

## License

Original analysis code in `code/` is licensed under the **MIT License**; see `LICENSE`.

Original derived tables and figures in `data/derived/`, `data/forensics/`, and `figures/` are licensed under **CC BY 4.0**; see `DATA_LICENSE.md`.

Underlying public datasets and third-party reference resources remain subject to their original terms of use.
