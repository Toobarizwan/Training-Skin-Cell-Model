#title: "preprocess_rds.rmd"

#This script converts a Seurat Object (.rds) into a pseudobulk object (this is a critical step to producing the individual-by-gene matrix which is necessary for the ranked-based percentiles).

## Did you activate the conda env?

## Load packages

library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(tibble)

## Read in .rds object

LD <- readRDS("/home/aliya/skin_cancer/MASLD_snRNA_seq_seurat_v4.rds")

## Filter cells and genes
# Keep cells with >200 genes
LD_filtered <- subset(LD, subset = nFeature_RNA > 200)

# Filter genes present in at least 3 cells
counts <- GetAssayData(LD_filtered, assay = "RNA", layer = "counts")
genesin3 <- rownames(counts)[Matrix::rowSums(counts > 0) >= 3]
LD_filtered <- LD_filtered[genesin3, ]

head(LD_filtered)

#Cell_type_detailed needs to be replaced with what your rds file is labelling the cell types
table(LD_filtered$Cell_type_detailed)

## Load Enformer metadata

enformer_metadata <- read.csv("/home/aliya/Liver/metadata.csv")
head(enformer_metadata)

# Get gene names
seurat_genes <- rownames(LD_filtered)
enformer_genes <- enformer_metadata$external_gene_name

# Find common genes
common_genes <- intersect(seurat_genes, enformer_genes)

# Debugging
cat("Genes in Seurat:", length(seurat_genes), "\n")
cat("Genes in Enformer:", length(enformer_genes), "\n")
cat("Common genes:", length(common_genes), "\n")
cat("Missing from Enformer:", length(setdiff(seurat_genes, enformer_genes)), "\n")
cat("Missing from Seurat:", length(setdiff(enformer_genes, seurat_genes)), "\n")

## Filter for common genes
counts <- GetAssayData(LD_filtered, assay = "RNA", layer = "counts")
common_genes <- intersect(rownames(counts), enformer_genes)
counts_sub <- counts[common_genes, , drop = FALSE]

## New object that is now filtered
LD_filteredEnf <- CreateSeuratObject(
  counts = counts_sub,
  meta.data = LD_filtered@meta.data)

## Split by cell type
head(LD_filteredEnf@meta.data)
table(LD_filteredEnf$Cell_type_detailed)

## Split into list by cell type

#Again, for your .rds file it will not be "Cell_type_detailed), so change it to be what the actual column name is

celltype_list <- SplitObject(LD_filteredEnf, split.by = "Cell_type_detailed")
cat("Cell types found:", length(celltype_list), "\n")
cat("Cell type names:\n")
print(names(celltype_list))

## Generate pseudobulk expression per cell type

for (ct in names(celltype_list)) {
  message("Processing cell type: ", ct)
  obj <- celltype_list[[ct]]

  # Check minimum cells
  if (ncol(obj) < 3) {
    message("Skipping ", ct, " - too few cells: ", ncol(obj))
    next
  }

  # Check DonorID column exists
  if (!"DonorID" %in% colnames(obj@meta.data)) {
    message("Skipping ", ct, " - DonorID column not found")
    next
  }

  # Normalize
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)

  # Aggregate by donor
  IbG <- AggregateExpression(
    obj,
    assays = "RNA",
    return.seurat = FALSE,
    group.by = "DonorID",
    verbose = FALSE
  )

  IbG_mat <- IbG$RNA

  # Convert to data frame with gene_name column
  IbG_df <- as.data.frame(IbG_mat)
  IbG_df <- tibble::rownames_to_column(IbG_df, var = "gene_name")

  # Clean cell type name for file path (remove spaces/special chars)
  ct_clean <- gsub("[^A-Za-z0-9_]", "_", ct)

  # Write TSV (Change to the directory you want it to be)
  out_path <- paste0("/home/tooba/Skin_Cell_Model/sample_file1/pseudobulk_", ct_clean, ".tsv")
  write.table(
    IbG_df,
    file = out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Written: ", out_path, " (", nrow(IbG_df), " genes x ", ncol(IbG_df)-1, " donors)")
}

#This is the pseudobulk data, it still needs to be preprocessed further for ctPred

