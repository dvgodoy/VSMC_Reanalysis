library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(reticulate)
library(scCustomize)
library(stringr)
library(readr)
library(clusterProfiler)
library(enrichplot)
library(DOSE)

set.seed(1234)
options(stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) User-configurable paths and parameters
# ─────────────────────────────────────────────────────────────────────────────
# In R Studio, retrieves the folder where this script is
# If you're not using R Studio, please set BASE_DIR manually
BASE_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(BASE_DIR)
META_FILE <- file.path(BASE_DIR, "metadata", "Carotid_MetaData_V1.txt")
RAW_DIR <- file.path(BASE_DIR, "expression", "raw")
MTX_FILE <- file.path(RAW_DIR, "Carotid_Expression_Matrix_raw_counts_V1.mtx.gz")
BARCODES_FILE <- file.path(RAW_DIR, "Carotid_Expression_Matrix_barcodes_V1.tsv.gz")
GENES_FILE <- file.path(RAW_DIR, "Carotid_Expression_Matrix_genes_V1.tsv.gz")
LATENT_FILE <- file.path(BASE_DIR, "embeddings", "latent.RData")

OUT_DIR <- file.path(BASE_DIR, "results")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
OBJ_DIR <- file.path(OUT_DIR, "objects")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

# Reticulate/scVI parameters.
RUN_SCVI <- FALSE
# If RUN_SCVI is FALSE, the arguments below are not used
# and it will use the latent.RData file from the embeddings folder instead
PYTHON_BIN <- "/path/to/python_env_with_scvi"
# Example
# PYTHON_BIN <- "/home/yourusername/miniconda3/envs/yourenv_with_scvi/bin/python3.12"
SCVI_BATCH_KEY <- "biosample_id"   # batch correction per experiment; biosample/library is closest here
SCVI_N_LATENT <- 50L
SCVI_MAX_EPOCHS <- 100L
UMAP_MIN_DIST <- 0.2
K_NEIGHBORS <- 15
LEIDEN_RESOLUTION <- 0.4

# ─────────────────────────────────────────────────────────────────────────────
# 2) Helper functions
# ─────────────────────────────────────────────────────────────────────────────
save_plot <- function(plot, filename, width = 9, height = 6, dpi = 300) {
  ggsave(file.path(FIG_DIR, filename), plot = plot, width = width, height = height, dpi = dpi)
}

sanitize_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

write_df <- function(x, filename) {
  write.csv(x, file.path(TAB_DIR, filename), row.names = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) Load counts and metadata safely
# ─────────────────────────────────────────────────────────────────────────────
print("Loading raw matrix and metadata")

raw_counts <- ReadMtx(
  mtx = MTX_FILE,
  cells = BARCODES_FILE,
  features = GENES_FILE,
  feature.column = 2,
  unique.features = TRUE
)

metadata <- read.delim(META_FILE, check.names = FALSE)
metadata <- metadata[metadata$NAME != "TYPE", , drop = FALSE]
rownames(metadata) <- metadata$NAME

# Convert numeric metadata columns from character to numeric.
numeric_cols <- c("n_umi", "n_genes", "percent_mito", "exon_prop", "entropy", "doublet_score")
for (cc in intersect(numeric_cols, colnames(metadata))) {
  metadata[[cc]] <- as.numeric(metadata[[cc]])
}

# Alignment is critical.
stopifnot(all(colnames(raw_counts) %in% rownames(metadata)))
metadata <- metadata[colnames(raw_counts), , drop = FALSE]
stopifnot(identical(colnames(raw_counts), rownames(metadata)))

metadata <- metadata %>%
  mutate(
    donor_id = factor(donor_id),
    biosample_id = factor(biosample_id),
    celltype = factor(celltype),
    condition = disease__ontology_label,
    condition2 = ifelse(condition == "normal", "control", "disease"),
    condition2 = factor(condition2, levels = c("control", "disease")),
    disease_subtype = case_when(
      condition == "normal" ~ "control",
      grepl("Symptomatic", biosample_id, ignore.case = TRUE) ~ "symptomatic",
      grepl("Asymptomatic", biosample_id, ignore.case = TRUE) ~ "asymptomatic",
      TRUE ~ "disease"
    ),
    disease_subtype = factor(disease_subtype, levels = c("control", "asymptomatic", "symptomatic", "disease"))
  )
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$NAME

write_df(as.data.frame(table(metadata$donor_id, metadata$condition2)), "donor_by_condition.csv")
write_df(as.data.frame(table(metadata$celltype, metadata$condition2)), "celltype_by_condition.csv")

# Keep all genes initially.
obj <- CreateSeuratObject(
  counts = raw_counts,
  project = "SCP2019_VSMC_Carotid",
  min.cells = 3,
  min.features = 200,
  meta.data = metadata
)
DefaultAssay(obj) <- "RNA"

print(paste0("Object contains ", ncol(obj), " nuclei and ", nrow(obj), " genes before gene-level filtering."))
print(table(obj$donor_id, obj$condition2))
print(table(obj$celltype, obj$condition2))

# ─────────────────────────────────────────────────────────────────────────────
# 4) QC and RNA-complexity mode diagnostics
# ─────────────────────────────────────────────────────────────────────────────
print("QC and complexity diagnostics")

# Metadata percent_mito is a fraction, so convert to percent for plots.
obj$percent.mt.metadata <- 100 * obj$percent_mito
obj$percent.mt.seurat <- PercentageFeatureSet(obj, pattern = "^MT-")
obj$percent.mt <- obj$percent.mt.metadata

obj$complexity <- obj$nFeature_RNA / obj$nCount_RNA
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

qc_fit <- lm(log10(nFeature_RNA) ~ log10(nCount_RNA), data = obj@meta.data)
obj$complexity_resid <- resid(qc_fit)

km <- kmeans(obj$complexity_resid, centers = 2, nstart = 50)
low_k <- names(which.min(tapply(obj$complexity_resid, km$cluster, mean)))
obj$qc_band <- ifelse(as.character(km$cluster) == low_k, "lower_complexity", "higher_complexity")
obj$qc_band <- factor(obj$qc_band, levels = c("higher_complexity", "lower_complexity"))

qc_summary <- obj@meta.data %>%
  group_by(donor_id, condition2, biosample_id, celltype, qc_band) %>%
  summarise(
    n = n(),
    median_nCount = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    median_percent_mt = median(percent.mt),
    median_complexity = median(complexity),
    median_complexity_resid = median(complexity_resid),
    .groups = "drop"
  )
write_df(qc_summary, "qc_summary_by_donor_celltype_band.csv")

qc_band_props <- obj@meta.data %>%
  count(celltype, qc_band) %>%
  group_by(celltype) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
write_df(qc_band_props, "qc_band_fraction_by_celltype.csv")

p_qc1 <- VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "complexity_resid"),
                 group.by = "condition2", pt.size = 0, ncol = 4)
save_plot(p_qc1, "qc_violin_by_condition.png", width = 14, height = 4)

p_qc2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "qc_band")
save_plot(p_qc2, "qc_nCount_vs_nFeature_qc_band.png", width = 7, height = 5)

p_qc3 <- ggplot(obj@meta.data, aes(x = complexity_resid, fill = qc_band)) +
  geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
  theme_classic() +
  labs(x = "Residual complexity: residuals from log10(nFeature) ~ log10(nCount)", y = "N nuclei")
save_plot(p_qc3, "qc_complexity_residual_histogram.png", width = 8, height = 5)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Normalization and variable genes
# ─────────────────────────────────────────────────────────────────────────────
print("Normalizing and selecting variable features")

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Exclude genes that are commonly technical/QC-driven from the HVG set.
technical_hvgs <- grep("^MT-|^RP[SL]|^HB[AB]|^MALAT1$", VariableFeatures(obj), value = TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj), technical_hvgs)
write_df(data.frame(variable_feature = VariableFeatures(obj)), "variable_features_after_technical_gene_removal.csv")
write_df(data.frame(removed_technical_hvg = technical_hvgs), "removed_technical_variable_features.csv")

top10 <- head(VariableFeatures(obj), 10)
p1_varf <- VariableFeaturePlot(obj)
p2_varf <- LabelPoints(plot = p1_varf, points = top10, repel = TRUE)
save_plot(p2_varf, "variable_features.png", width = 8, height = 5)

# ─────────────────────────────────────────────────────────────────────────────
# 6) scVI embedding + Seurat graph/UMAP/clustering
# ─────────────────────────────────────────────────────────────────────────────
print("Dimensionality reduction and clustering")

if (RUN_SCVI) {
  Sys.setenv(RETICULATE_PYTHON = PYTHON_BIN)
  use_python(PYTHON_BIN, required = TRUE)

  scvi <- import("scvi", convert = FALSE)
  scipy <- import("scipy", convert = FALSE)

  adata <- as.anndata(
    obj,
    file_path = OBJ_DIR,
    file_name = "carotid_scvi_input.h5ad",
    assay = "RNA",
    main_layer = "counts",
    other_layers = "data"
  )

  obs_cols <- py_to_r(adata$obs$columns$tolist())

  # Ensure X is a sparse matrix. as.anndata(..., main_layer "counts") should place raw counts in X.
  adata$X <- scipy$sparse$csr_matrix(adata$X)
  
  scvi$settings$dl_num_workers <- 0L
  scvi$settings$seed <- 1234L

  scvi$model$SCVI$setup_anndata(adata, batch_key = SCVI_BATCH_KEY)
  model <- scvi$model$SCVI(adata, n_latent = SCVI_N_LATENT)
  model$train(max_epochs = SCVI_MAX_EPOCHS, early_stopping = TRUE)

  latent <- py_to_r(model$get_latent_representation())
  latent <- as.matrix(latent)
  rownames(latent) <- colnames(obj)
  colnames(latent) <- paste0("scvi_", seq_len(ncol(latent)))

  save(latent, file=LATENT_FILE)
  save(model, file='model_scvi_50.RData')  
} else {
  load(LATENT_FILE)
}

obj[["scvi"]] <- CreateDimReducObject(embeddings = latent, assay = "RNA")
dims_use <- seq_len(ncol(latent))

obj <- FindNeighbors(obj, reduction = "scvi", dims = dims_use, k.param = K_NEIGHBORS, verbose = FALSE)
obj <- FindClusters(obj, algorithm = 4, resolution = LEIDEN_RESOLUTION, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "scvi", dims = dims_use, min.dist = UMAP_MIN_DIST, verbose = FALSE)

# Majority-vote cell-type label for each Leiden cluster, while preserving original metadata celltype.
cluster_to_celltype <- obj@meta.data %>%
  count(seurat_clusters, celltype, name = "n") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(as.character(seurat_clusters)))

write_df(cluster_to_celltype, "cluster_majority_celltype.csv")

cluster_map <- setNames(as.character(cluster_to_celltype$celltype), as.character(cluster_to_celltype$seurat_clusters))
tmp <- cluster_map[as.character(obj$seurat_clusters)]
names(tmp) <- colnames(obj)
obj$cluster_majority_celltype <- tmp

p_umap_celltype <- DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE, pt.size = 0.35) + NoLegend()
p_umap_donor <- DimPlot(obj, reduction = "umap", group.by = "donor_id", pt.size = 0.35)
p_umap_condition <- DimPlot(obj, reduction = "umap", group.by = "condition2", pt.size = 0.35)
p_umap_qc <- DimPlot(obj, reduction = "umap", group.by = "qc_band", pt.size = 0.35)
save_plot((p_umap_celltype | p_umap_donor) / (p_umap_condition | p_umap_qc), "umap_overview.png", width = 13, height = 10)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Marker genes and paper-focused visualization
# ─────────────────────────────────────────────────────────────────────────────
print("Marker genes and paper-focused plots")

Idents(obj) <- "celltype"
celltype_markers_sc <- FindAllMarkers(
  obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)
top_markers = celltype_markers_sc %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
write_df(top_markers %>% as.data.frame(), 'top_markers.csv')
write_df(celltype_markers_sc, "single_cell_celltype_markers_exploratory.csv")

paper_genes <- c(
  "HDAC9", "CD68", "LGALS3", "SERPINE2", "MALAT1",
  "ACTA2", "TAGLN", "MYH11", "CNN1", "SMTN",
  "PDE4D", "ANK3", "LAMA2",
  "PECAM1", "VWF", "CDH5",
  "DCN", "LUM", "COL1A1", "COL6A3",
  "PTPRC", "CD3D", "CD3E", "IL7R", "LYZ", "MSR1", "MRC1"
)
paper_genes <- intersect(paper_genes, rownames(obj))
write_df(data.frame(gene = paper_genes), "paper_genes_present.csv")

p_dot <- DotPlot(obj, features = paper_genes, group.by = "celltype", split.by = "condition2", cols = c("red", "blue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8))
save_plot(p_dot, "dotplot_paper_genes_by_celltype_condition.png", width = 15, height = 6)

p_dot2 <- DotPlot(obj, features = top_markers$gene, group.by = "celltype", split.by = "condition2", cols = c("red", "blue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8))
save_plot(p_dot2, "dotplot_top_markers_by_celltype_condition.png", width = 15, height = 6)

p_feat <- FeaturePlot(obj, features = intersect(c("HDAC9", "CD68", "LGALS3", "SERPINE2", "ACTA2", "MYH11"), rownames(obj)), ncol = 3)
save_plot(p_feat, "featureplot_key_genes.png", width = 12, height = 8)

# Donor-level summaries for target genes.
target_expr <- lapply(paper_genes, function(g) {
  x <- FetchData(obj, vars = c(g, "donor_id", "biosample_id", "condition2", "celltype"))
  colnames(x)[1] <- "expr"
  x %>%
    group_by(gene = g, donor_id, biosample_id, condition2, celltype) %>%
    summarise(
      n_cells = n(),
      mean_log_norm_expr = mean(expr),
      pct_expr = mean(expr > 0),
      .groups = "drop"
    )
}) %>% bind_rows()
write_df(target_expr, "target_gene_expression_by_donor_celltype.csv")

# ─────────────────────────────────────────────────────────────────────────────
# 8) Cell-type composition summaries
# ─────────────────────────────────────────────────────────────────────────────
print("Cell-type composition summaries")

celltype_counts <- obj@meta.data %>%
  count(donor_id, condition2, disease_subtype, celltype, name = "n_cells") %>%
  group_by(donor_id) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()
write_df(celltype_counts, "celltype_composition_by_donor_biosample.csv")

p_comp <- ggplot(celltype_counts, aes(x = donor_id, y = prop, fill = celltype)) +
  geom_col() +
  facet_grid(. ~ condition2, scales = "free_x", space = "free_x") +
  theme_classic() +
  labs(x = "Donor", y = "Fraction of nuclei", fill = "Cell type")
save_plot(p_comp, "celltype_composition_by_donor.png", width = 9, height = 5)

# ─────────────────────────────────────────────────────────────────────────────
# 9) Exploratory single-cell disease-vs-control tests
# ─────────────────────────────────────────────────────────────────────────────
print("Exploratory single-cell disease-vs-control DE")

run_singlecell_de <- function(obj, celltype_name) {
  obj_ct <- subset(obj, subset = celltype == celltype_name)
  
  Idents(obj_ct) <- "condition2"
  latent_vars <- intersect(c("nCount_RNA", "percent.mt", "complexity_resid"), colnames(obj_ct@meta.data))

  # avg_log2FC > 0  # higher in disease
  # avg_log2FC < 0  # higher in control
  res <- FindMarkers(
    obj_ct,
    ident.1 = "disease",
    ident.2 = "control",
    assay = "RNA",
    test.use = "LR",
    latent.vars = latent_vars,
    min.pct = 0.05,
    logfc.threshold = 0,
    verbose = FALSE
  ) %>%
    rownames_to_column("gene") %>%
    mutate(
      celltype = celltype_name,
      contrast = "disease_vs_control",
      method = "Seurat_LR_cell_level_exploratory",
    ) %>%
    relocate(celltype, gene, contrast)

  write_df(res, paste0("singlecell_DE_disease_vs_control_", sanitize_id(celltype_name), ".csv"))
  res
}

celltypes_to_test = c("Endothelial", "Lymphocyte", "Macrophage", "Pericyte", "Vascular Smooth Muscle")

# Takes quite a while to run!
sc_de_results <- lapply(celltypes_to_test, function(ct) {
  run_singlecell_de(obj, ct)
}) %>% bind_rows()

write_df(sc_de_results, "singlecell_DE_disease_vs_control_all_celltypes_exploratory.csv")

# ─────────────────────────────────────────────────────────────────────────────
# 10) Gene Set Enrichment
# ─────────────────────────────────────────────────────────────────────────────
print("Gene set enrichment")

for (ctype in unique(celltype_markers_sc$cluster)) {
  tmp_df <- subset(celltype_markers_sc, cluster == ctype)

  original_gene_list <- tmp_df$avg_log2FC
  names(original_gene_list) <- rownames(tmp_df)
  names(original_gene_list) <- tmp_df$gene
  gene_list <- na.omit(original_gene_list)
  gene_list <- sort(gene_list, decreasing=T)

  organism <- 'org.Hs.eg.db'
  gse <- gseGO(geneList=gene_list, ont='ALL', keyType='SYMBOL', nPerm=10000, minGSSize = 3, maxGSSize = 800, 
              pvalueCutoff = 0.05, verbose = T, OrgDb=organism, pAdjustMethod = 'none')

  act_sup_plot <- dotplot(gse, showCategory=10, split='.sign') + facet_grid(.~.sign)
  save_plot(act_sup_plot, paste0('act_sup_', ctype, '.png'), width=10, height=8)

  x <- pairwise_termsim(gse)
  netw_plot <- emapplot(x, showCategory=10)
  save_plot(netw_plot, paste0('netw_', ctype, '.png'), width=10, height=8)

  ridge_cat_plot <- ridgeplot(gse, showCategory=10) + labs(x='Enrichment Distribution')
  save_plot(ridge_cat_plot, paste0('ridge_', ctype, '.png'), width=10, height=8)
}

# ─────────────────────────────────────────────────────────────────────────────
# 11) Additional plots
# ─────────────────────────────────────────────────────────────────────────────

keep_features <- function(object, genes) {
  genes <- unique(genes)
  genes[genes %in% rownames(object)]
}

keep_feature_list <- function(object, gene_list) {
  gene_list <- lapply(gene_list, function(x) keep_features(object, x))
  gene_list[lengths(gene_list) > 0]
}

# Celltype-condition group for clean dotplots
cell_order <- c(
  "Endothelial",
  "Fibroblast",
  "Lymphocyte",
  "Macrophage",
  "Pericyte",
  "Unknown1",
  "Unknown2",
  "Vascular Smooth Muscle"
)

condition_order <- c("control", "disease")

wanted_levels <- as.vector(
  outer(cell_order, condition_order, paste, sep = " | ")
)

obj$celltype_condition <- paste(obj$celltype, obj$condition2, sep = " | ")
obj$celltype_condition <- factor(
  obj$celltype_condition,
  levels = wanted_levels[wanted_levels %in% unique(obj$celltype_condition)]
)

paper_gene_sets <- list(
  # https://www.nature.com/articles/s41467-018-03394-7
  # SMARCA4 = BRG1, SERPINE2 = PN1
  HDAC9_axis = c(
    "HDAC9", "MALAT1", "SMARCA4", "SERPINE2"
  ),
  # https://www.ahajournals.org/doi/10.1161/ATVBAHA.121.316600
  # TAGLN = SM22
  # https://www.babraham.ac.uk/sites/default/files/media/files/30385745.pdf
  VSMC_contractile = c(
    "ACTA2", "TAGLN", "MYH11", "CNN1", "SMTN"
  ),
  # Macrophage-Like vSMCs
  VSMC_immune_like = c(
    "CD68", "LGALS3", "LYZ", "MSR1", "APOE"
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9533272/
  Endothelial_controls = c(
    "PECAM1", "VWF", "CDH5", "CLDN5", "KDR", "FLT1", "TEK"
  )
)

paper_gene_sets <- keep_feature_list(obj, paper_gene_sets)

p_paper <- DotPlot(
  obj,
  features = paper_gene_sets,
  group.by = "celltype_condition",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle("HDAC9 / VSMC phenotype panel")

save_plot(p_paper, "dotplot_paper_focused_HDAC9_VSMC_panel.png", width = 15, height = 8)

vsmc <- subset(obj, subset = celltype == "Vascular Smooth Muscle")

vsmc$donor_condition <- paste(vsmc$donor_id, vsmc$condition2, sep = " | ")

vsmc_gene_sets <- list(
  Paper_axis = c(
    "HDAC9", "MALAT1", "SMARCA4", "SERPINE2", "CD68", "LGALS3"
  ),
  Contractile = c(
    "ACTA2", "TAGLN", "MYH11", "CNN1", "SMTN"
  ),
  Disease_higher_exploratory = vsmc_disease_higher_exploratory,
  Control_higher_exploratory = vsmc_control_higher_exploratory
)

vsmc_gene_sets <- keep_feature_list(vsmc, vsmc_gene_sets)

p_vsmc_de <- DotPlot(
  vsmc,
  features = vsmc_gene_sets,
  group.by = "donor_condition",
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("VSMC donor-level disease/control candidate genes")

save_plot(p_vsmc_de, "dotplot_VSMC_disease_control_candidates_by_donor.png", width = 14, height = 5)

annotation_sets <- list(
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "KDR", "FLT1", "TEK"),
  VSMC = c("ACTA2", "TAGLN", "MYH11", "CNN1", "SMTN"),
  Fibroblast = c("DCN", "LUM", "COL1A1", "COL1A2", "COL6A3", "MFAP5"),
  Pericyte = c("RGS5", "ABCC9", "PDGFRB", "CSPG4", "MCAM"),
  Macrophage = c("LYZ", "CD68", "CD163", "MSR1", "MRC1", "C1QA", "C1QB", "C1QC"),
  Lymphocyte = c("PTPRC", "CD2", "CD3D", "CD3E", "TRAC", "IL7R", "BCL11B")
)

annotation_sets <- keep_feature_list(obj, annotation_sets)

p_annotation <- DotPlot(
  obj,
  features = annotation_sets,
  group.by = "celltype",
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("Cell-type annotation marker sanity check")

save_plot(p_annotation, "dotplot_celltype_annotation_sanity_check.png", width = 14, height = 6)

aggregate.vsmc <- AggregateExpression(obj, assays = "RNA", group.by = c("celltype", "condition2"), return.seurat = TRUE, verbose = FALSE)

Idents(obj) = paste0(obj$celltype, sep='_', obj$condition2)
vsmc_markers <- FindMarkers(obj, ident.1 = "Vascular Smooth Muscle_control", ident.2 = "Vascular Smooth Muscle_disease", assay = "RNA", slot = "data", verbose = FALSE)
endot_markers <- FindMarkers(obj, ident.1 = "Endothelial_disease", ident.2 = "Endothelial_control", assay = "RNA", slot = "data", verbose = FALSE)

genes.to.label <- rownames(vsmc_markers)[1:15]
p1 <- CellScatter(aggregate.vsmc, cell1 = "Vascular Smooth Muscle_control", cell2 = "Vascular Smooth Muscle_disease", highlight = genes.to.label) +
  ggtitle("Vascular Smooth Muscle Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
genes.to.label <- rownames(endot_markers)[1:15]
p3 <- CellScatter(aggregate.vsmc, cell1 = "Endothelial_control", cell2 = "Endothelial_disease", highlight = genes.to.label) +
  ggtitle("Endothelial Cells")
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
save_plot(p1 + p3, "genes_across_conditions.png", width = 12, height = 6)

genes.to.label = paper_genes
p1 <- CellScatter(aggregate.vsmc, cell1 = "Vascular Smooth Muscle_control", cell2 = "Vascular Smooth Muscle_disease", highlight = genes.to.label) +
  ggtitle("Vascular Smooth Muscle Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p3 <- CellScatter(aggregate.vsmc, cell1 = "Endothelial_control", cell2 = "Endothelial_disease", highlight = genes.to.label) +
  ggtitle("Endothelial Cells")
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
save_plot(p1 + p3, "paper_genes_across_conditions.png", width = 12, height = 6)

Idents(obj) = obj$celltype

genes <- c("CD68", "LGALS3", "HDAC9")

celltypes_keep <- setdiff(
  sort(unique(obj$celltype)),
  c("Unknown1", "Unknown2", "Fibroblast")
)

make_dot <- function(seu, title) {
  seu$celltype <- factor(seu$celltype, levels = celltypes_keep)
  
  DotPlot(
    seu,
    features = genes,
    group.by = "celltype",
    cols = c("grey90", "blue"),
    col.min = -1,
    col.max = 2.5,
    dot.min = 0,
    scale.min = 0,
    scale.max = 100
  ) +
    ggtitle(title) +
    RotatedAxis() +
    guides(
      colour = guide_colorbar(
        title = "Avg. expr.",
        title.position = "top",
        barwidth = unit(2.2, "cm"),
        barheight = unit(0.25, "cm")
      ),
      size = guide_legend(
        title = "% expr.",
        title.position = "top",
        nrow = 1,
        override.aes = list(colour = "grey40")
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.25, "cm"),
      legend.spacing.x = unit(0.15, "cm"),
      legend.margin = margin(t = -2, r = 0, b = 0, l = 0),
    )
}

obj_control <- subset(
  obj,
  subset = !(celltype %in% c("Unknown1", "Unknown2", "Fibroblast")) &
    condition2 == "control"
)

obj_disease <- subset(
  obj,
  subset = !(celltype %in% c("Unknown1", "Unknown2", "Fibroblast")) &
    condition2 == "disease"
)

p1 <- make_dot(obj_control, "Control")
p2 <- make_dot(obj_disease, "Disease")
combined <- p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(combined, 'three_genes.png', width=12, height=6)

# ─────────────────────────────────────────────────────────────────────────────
# 12) Final save and session info
# ─────────────────────────────────────────────────────────────────────────────
print("Saving outputs")

saveRDS(obj, file.path(OBJ_DIR, "carotid_seurat_improved.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

print(paste0("Done. Main outputs written to: ", OUT_DIR))