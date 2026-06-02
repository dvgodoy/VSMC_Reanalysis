# VSMC Reanalysis

Reanalysis of the single-nucleus RNA-seq dataset from **“Vascular smooth muscle cell phenotype switching in carotid atherosclerosis”** by Chou et al. This project focuses on reproducing and extending the paper’s main biological observation: vascular smooth muscle cells (VSMCs) in carotid atherosclerosis show evidence of phenotype modulation, including altered contractile identity and expression of immune/macrophage-associated genes.

The analysis is implemented in R, with optional use of Python/scVI through `reticulate` for batch-aware dimensionality reduction.

---

## 1. Project goals

This repository reanalyzes the raw count data associated with the original study. The main goals are:

1. Reconstruct a Seurat object from the raw 10x-style expression matrix and metadata.
2. Perform quality-control diagnostics, including RNA-complexity mode estimation.
3. Use the provided scVI latent representation, or optionally regenerate it, to obtain a batch-aware embedding.
4. Reproduce the major cell-type structure reported in the paper.
5. Focus on VSMC phenotype switching, especially the HDAC9-associated mechanism discussed by the authors.
6. Examine expression of contractile VSMC genes, immune/macrophage-like markers, and genes involved in the proposed HDAC9–MALAT1–BRG1/SERPINE2 axis.
7. Perform exploratory disease-versus-control comparisons while explicitly accounting for the limitations caused by donor/condition confounding.
8. Run gene set enrichment analysis to connect VSMC marker genes with biological-process ontologies, especially muscle cell differentiation and contractile programs.

This project is intended as a transparent reanalysis and educational reproduction, not as a replacement for the original publication.

---

## 2. Data sources

### Original data

The original single-nucleus RNA-seq data were obtained from the Broad Institute Single Cell Portal study page:

- [Single Cell Portal SCP2019](https://singlecell.broadinstitute.org/single_cell/study/SCP2019/)

Only the **raw expression data** were used in this reanalysis. The normalized matrix provided by the original authors was not used, because this project reconstructs normalization, quality control, embedding, clustering, marker analysis, and visualization from raw counts.

### Copy of raw data in this repository

For convenience and reproducibility, a copy of the raw expression files used here is included in:

- [`expression/raw`](expression/raw)

Expected raw input files:

```text
expression/raw/
├── Carotid_Expression_Matrix_raw_counts_V1.mtx.gz
├── Carotid_Expression_Matrix_barcodes_V1.tsv.gz
└── Carotid_Expression_Matrix_genes_V1.tsv.gz
```

The cell-level metadata are stored in:

```text
metadata/
└── Carotid_MetaData_V1.txt
```

The repository also includes a precomputed scVI latent embedding:

```text
embeddings/
└── latent.RData
```

By default, the analysis uses this file instead of retraining scVI.

---

## 3. Repository structure

```text
VSMC_Reanalysis/
├── embeddings/
│   └── latent.RData
├── expression/
│   └── raw/
│       ├── Carotid_Expression_Matrix_raw_counts_V1.mtx.gz
│       ├── Carotid_Expression_Matrix_barcodes_V1.tsv.gz
│       └── Carotid_Expression_Matrix_genes_V1.tsv.gz
├── metadata/
│   └── Carotid_MetaData_V1.txt
├── results/
│   ├── figures/
│   ├── objects/
│   ├── tables/
│   └── sessionInfo.txt
├── project.R
├── requirements.txt
├── file_supplemental_info.tsv
├── Reanalysis of Vascular Smooth Muscle Cell Phenotype Switching in Carotid Atherosclerosis.pdf
├── LICENSE
└── README.md
```

The main analysis is contained in:

```text
project.R
```

---

## 4. Software requirements

### R

The analysis was written in R and uses the following packages:

```r
Seurat
Matrix
dplyr
tidyr
tibble
ggplot2
patchwork
reticulate
scCustomize
stringr
readr
clusterProfiler
enrichplot
DOSE
org.Hs.eg.db
rstudioapi
```

A typical installation command is:

```r
install.packages(c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "reticulate",
  "stringr",
  "readr",
  "rstudioapi"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "org.Hs.eg.db"
))
```

### Python (only needed to regenerate scVI embeddings)

The repository includes a precomputed scVI latent representation. Therefore, Python is not required for the default run if:

```r
RUN_SCVI <- FALSE
```

To regenerate the scVI embedding, install a Python environment containing:

```text
scvi-tools
scanpy
scipy
anndata
numpy
pandas
```

Example using conda or mamba:

```bash
conda create -n scvi_env python=3.12 -y
conda activate scvi_env
pip install scvi-tools scanpy scipy anndata numpy pandas
```

Then update the `PYTHON_BIN` variable in `project.R`:

```r
PYTHON_BIN <- "/path/to/scvi_env/bin/python"
RUN_SCVI <- TRUE
```

---

## 5. Running the analysis

Clone the repository:

```bash
git clone https://github.com/dvgodoy/VSMC_Reanalysis.git
cd VSMC_Reanalysis
```

Open `project.R` in RStudio and run the script.

The script uses:

```r
BASE_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(BASE_DIR)
```

If running outside RStudio, manually set `BASE_DIR` near the top of `project.R`:

```r
BASE_DIR <- "/path/to/VSMC_Reanalysis"
setwd(BASE_DIR)
```

Then run:

```bash
Rscript project.R
```

The default configuration is:

```r
RUN_SCVI <- FALSE
```

This loads the precomputed latent embedding from:

```text
embeddings/latent.RData
```

### Reproducing the scVI embedding

To regenerate the embedding:

1. Install a Python environment with `scvi-tools`.
2. Set `PYTHON_BIN` in `project.R`.
3. Set:

```r
RUN_SCVI <- TRUE
```

4. Run the full script.

The scVI configuration used in the script is:

```r
SCVI_BATCH_KEY <- "biosample_id"
SCVI_N_LATENT <- 50L
SCVI_MAX_EPOCHS <- 100L
```

The embedding is then used for:

```r
FindNeighbors(reduction = "scvi")
FindClusters(algorithm = 4, resolution = 0.4)
RunUMAP(reduction = "scvi", min.dist = 0.2)
```

---

## 6. Main outputs

### Figures

The main figures are written to:

```text
results/figures/
```

Selected important figures:

```text
qc_violin_by_condition.png
qc_nCount_vs_nFeature_qc_band.png
qc_complexity_residual_histogram.png
umap_overview.png
celltype_composition_by_donor.png
dotplot_celltype_annotation_sanity_check.png
dotplot_paper_genes_by_celltype_condition.png
dotplot_paper_focused_HDAC9_VSMC_panel.png
dotplot_VSMC_disease_control_candidates_by_donor.png
featureplot_key_genes.png
act_sup_Vascular Smooth Muscle.png
ridge_Vascular Smooth Muscle.png
```

### Tables

The main tables are written to:

```text
results/tables/
```

Selected important tables:

```text
donor_by_condition.csv
celltype_by_condition.csv
qc_summary_by_donor_celltype_band.csv
qc_band_fraction_by_celltype.csv
variable_features_after_technical_gene_removal.csv
removed_technical_variable_features.csv
single_cell_celltype_markers_exploratory.csv
top_markers.csv
cluster_majority_celltype.csv
celltype_composition_by_donor_biosample.csv
target_gene_expression_by_donor_celltype.csv
singlecell_DE_disease_vs_control_all_celltypes_exploratory.csv
singlecell_DE_disease_vs_control_Vascular_Smooth_Muscle.csv
```

### Seurat object

The processed Seurat object is saved to:

```text
results/objects/carotid_seurat_improved.rds
```

### Session information

R session information is saved to:

```text
results/sessionInfo.txt
```

---

## 7. Citation

If you use this repository, cite the original paper:

> Chou EL, Lino Cardenas CL, Chaffin M, Arduini AD, Juric D, Stone JR, LaMuraglia GM, Eagleton MJ, Conrad MF, Isselbacher EM, Ellinor PT, Lindsay ME. **Vascular smooth muscle cell phenotype switching in carotid atherosclerosis.** *JVS: Vascular Science*. 2022;3:41–47. https://doi.org/10.1016/j.jvssci.2021.11.002

Also cite the original data source:

- [Single Cell Portal SCP2019](https://singlecell.broadinstitute.org/single_cell/study/SCP2019/)
