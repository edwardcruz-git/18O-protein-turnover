# Proteome-wide quantification of protein turnover in frog and fly embryos reveals divergent strategies of maternal inheritance

Reproducible analysis supplement for the paper of the same title.

This repository contains the data, intermediate files, reference tables, and analysis scripts underlying the paper. Each analysis is provided both as a rendered document, so the code can be read alongside the output it produced, and as its source `.Rmd`, so the analysis can be re-run.

## Data availability

The mass spectrometry proteomics data have been deposited to the ProteomeXchange consortium via the PRIDE partner repository under accessions **PXD081127** (DDA) and **PXD084487** (DIA). Compressed copies of the raw search outputs are also included here under `raw/Data/` so the full pipeline can be re-run without retrieving the deposit separately.

**This is a large repository (approximately 800 MB).** It bundles the compressed raw mass spectrometry files and DIA-NN reports so the analysis is self-contained, which makes cloning slower than a code-only repository. If you only want to read the code and results, browse the rendered documents in `analysis/` directly on GitHub without cloning.

## Repository structure

- `analysis/` — rendered analysis documents (`.md`) with their figures. This is what to read to follow the analysis.
- `scripts/` — the source `.Rmd` files. Use these to re-run the analysis.
- `raw/` — all inputs the scripts read. Contains `Data/` (raw search outputs and derived intermediates), `Files/` (fits, reference tables, absolute concentrations, FASTA databases), and `Systems/` (InterProScan and IUPred2A annotation).

## Analyses

The documents are listed in suggested execution order. Because every intermediate file is included in `raw/`, each document finds its inputs already present and can also be run on its own, so the order below is a guide to the pipeline rather than a strict requirement. See the notes at the end for one exception where two documents depend on each other.

| # | Analysis | Description | Rendered |
|---|----------|-------------|----------|
| 1 | Frog_GB_NYS-Proteomics | *Xenopus* normalization set and yolk-normalized dataset | [view](analysis/Frog_GB_NYS-Proteomics.md) |
| 2 | Frog_AA_Incorporation_Early-T1 | *Xenopus* early (2-cell) timeseries — amino acid incorporation and raw peptide filtering | [view](analysis/Frog_AA_Incorporation_Early-T1.md) |
| 3 | Frog_AA_Incorporation_MBT-T2 | *Xenopus* gastrulation timeseries — amino acid incorporation and raw peptide filtering | [view](analysis/Frog_AA_Incorporation_MBT-T2.md) |
| 4 | Fly_AA_Incorporation_MBT-T2 | *Drosophila* gastrulation timeseries — amino acid incorporation and raw peptide filtering | [view](analysis/Fly_AA_Incorporation_MBT-T2.md) |
| 5 | Frog_2C_SimpModel | *Xenopus* 2-cell — protein decay model fitting | [view](analysis/Frog_2C_SimpModel.md) |
| 6 | Frog_Gast_SimpModel | *Xenopus* gastrulation — protein decay model fitting | [view](analysis/Frog_Gast_SimpModel.md) |
| 7 | Fly_Gast_SimpModel | *Drosophila* gastrulation — protein decay model fitting | [view](analysis/Fly_Gast_SimpModel.md) |
| 8 | turnover_abs_quant | Absolute protein concentrations for frog and fly | [view](analysis/turnover_abs_quant.md) |
| 9 | Comp_Frog_Early-MBT | Comparison of *Xenopus* 2-cell and gastrulation turnover | [view](analysis/Comp_Frog_Early-MBT.md) |
| 10 | AA-Tracking_Fly_and_Frog | Amino acid supply and demand across frog and fly | [view](analysis/AA-Tracking_Fly_and_Frog.md) |
| 11 | Turnover_Species_Comparison | Frog versus fly turnover across orthologues | [view](analysis/Turnover_Species_Comparison.md) |
| 12 | fly_turnover_analysis | *Drosophila* turnover with structural and functional annotation | [view](analysis/fly_turnover_analysis.md) |
| 13 | Frog_Early_Analysis | *Xenopus* turnover with structural and functional annotation | [view](analysis/Frog_Early_Analysis.md) |
| 14 | turnover_and_absolute | Relationship between turnover and absolute abundance | [view](analysis/turnover_and_absolute.md) |

## Reproducing the analysis

The scripts read their inputs using paths relative to `raw/` (for example `Data/XLA_O18/...` and `Files/Reference/...`). Two things follow from this and both are required, or the scripts will not find their inputs.

**1. Run with the working directory set to `raw/`.** The `.Rmd` files live in `scripts/`, but their paths are written relative to `raw/`. Because R Markdown defaults the working directory to the file's own location, you must set the knit root directory to `raw/` before rendering. Set it in a document with

```r
knitr::opts_knit$set(root.dir = "path/to/raw")
```

For running interactively rather than knitting, `setwd()` to `raw/` before sourcing.

**2. Decompress the raw files first.** The large raw search outputs are stored compressed as `.xz` under `raw/Data/`, but the scripts read them by their `.csv` names. Decompress them before running so the expected files exist. From inside `raw/`

```bash
find . -name "*.xz" -exec unxz {} \;
```

`unxz` strips the `.xz` extension and leaves the `.csv` the scripts expect. Only the documents that read raw data need this (the amino acid incorporation documents and the GB proteomics document). The downstream documents run from the derived intermediates already included uncompressed.

Session and package versions are recorded at the end of each rendered document.

## Notes

**Run order and the GB / Early-T1 dependency.** The normalization approach was finalized late in the project and is validated against the metabolomic labeling data, so `Frog_GB_NYS-Proteomics` and `Frog_AA_Incorporation_Early-T1` read from each other's outputs and are mutually dependent rather than strictly sequential. Because all intermediate files are included in `raw/`, this does not affect re-running, each document finds the inputs it needs already present. To regenerate both from scratch, run `Frog_AA_Incorporation_Early-T1` up to the point that produces the amino acid matrix, then run `Frog_GB_NYS-Proteomics`, then complete the remaining steps.

**"MBT" in script and file names.** Some scripts, figure folders, and intermediate files are labeled `MBT` for historical reasons. An earlier version of the study measured a different developmental window, and although the manuscript text uses the updated naming (gastrulation), the scripts were never renamed. Any `MBT` label in this repository refers to the gastrulation timeseries.

**FASTA files.** Protein search databases are split across `raw/Files/Reference/` and `raw/Files/FASTA/` (with a `DIA-NN/` subfolder) for historical reasons. The frog database was moved to `Files/FASTA/` when it was updated mid-project and the fly database was left in `Files/Reference/`, so scripts read from both. The current frog and fly databases are included in both folders for completeness, and the contaminant-containing versions under `DIA-NN/` are those used for the DIA-NN searches, which require them. Database choice does not affect the reported results.
