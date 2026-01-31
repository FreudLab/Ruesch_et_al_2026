# IL-15 Signaling During a Critical NK Cell Developmental Window Subverts Maturation and KIR Acquisition

![R](https://img.shields.io/badge/R-4.4.0-blue.svg)
![Python](https://img.shields.io/badge/Python-3.12.7-yellow.svg)

## Overview
This repository contains the code and analysis scripts for the study **"IL-15 Signaling During a Critical NK Cell Developmental Window Subverts Maturation and KIR Acquisition"** (Ruesch et al., 2026).

This study investigates the transcriptional and epigenetic trajectories of developing human NK cells under distinct IL-15 cytokine signaling schedules. The repository includes:
1.  **CITE-seq Analysis (Python):** 5' Capture single-cell analysis of transcript and protein expression from FACS-purified Innate Lymphoid Cells (ILCs).
2.  **Methylation Analysis (R):** Processing and analysis of methylation array data from FACS-sorted NK subsets.

## Data Availability
* **CITE-seq Data:** Raw and processed files are available at GEO under accession [GSE318250](https://www.ncbi.nlm.nih.gov/geo/).
* **Methylation Data:** Raw .idat and processed matrix files are available at GEO under accession [GSE318251](https://www.ncbi.nlm.nih.gov/geo/).

## Folder Structure
* `Processing/`: Scripts for initial processing of raw .idat files (R) and generating CITE-seq count matrices (Python).
* `Analysis/`: 
    * `CITEseq_clustering.ipynb`: Dimensionality reduction, clustering, and annotation.
    * `Methylation_RnBeads.R`: Differential methylation analysis.
* `Figures/`: Code to reproduce specific figure panels (e.g., Figure 2, Figure S5).
* `Processed_Data/`: Completely processed mudata and RnBeadSets used in the article
* `Envs/`: Contains the `environment.yml` and `sessionInfo.txt` files for reproducibility.

## Computational Requirements

### 1. Python (CITE-seq Analysis)
Analysis was performed using **Python v3.12.7**. 
Key libraries used include:
* **Scanpy** v1.11.5
* **Muon** v0.1.6
* **Decoupler** v1.8.0
* **Leidenalg** v0.10.2
* **Matplotlib** v3.10.0
* **Seaborn** v0.13.2

*To reproduce the exact environment, please refer to `Envs/environment.yml`.*

### 2. R (Methylation Analysis)
Analysis was performed using **R v4.4.0**.
Key libraries used include:
* **RnBeads** (BioConductor)
* **ggplot2**

**Note on Visualization:** While statistical processing of methylation data was performed in R using RnBeads, final figure aesthetic adjustments and plots were generated using **GraphPad Prism** using the exported processed tables found in the `Output/` folder.

## Citation

If you used this data or code, please cite this work as:

> **Ruesch, M.A., et al. (2026).** IL-15 Signaling During a Critical NK Cell Developmental Window Subverts Maturation and KIR Acquisition. *[Under Review]*.

If you use the code provided in this repository, please also reference the repository URL:
> https://github.com/FreudLab/NK-Methylation-Development

## Code Environments
To replicate the Python environment used for the CITE-seq analysis:

```bash
conda env create -f Envs/environment.yml
conda activate NK_Dev_2026
```

To reproduce the exact R environment used for the methylation analysis:

1. Install [renv](https://rstudio.github.io/renv/articles/renv.html).
2. Download the `renv.lock` file from this repository.
3. Run the following in R:
   ```R
   renv::restore()
