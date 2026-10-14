#we want to find correlations betwwn cell types and create a heatmap

library(dplyr)
library(readr)

pseudo_dir <- "/home/tooba/Skin_Cell_Model/pseudobulk/"

# Start with just Arterial_EC to confirm this works before looping
df <- read_tsv(file.path(pseudo_dir, "pseudobulk_Arterial_EC.tsv"))

# All columns except gene_name are donor sample columns
donor_cols <- setdiff(colnames(df), "gene_name")

# Mean expression per gene, across all donors, for this cell type
arterial_ec_mean <- df %>%
  rowwise() %>%
  mutate(mean_expression = mean(c_across(all_of(donor_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(gene_name, mean_expression)

head(arterial_ec_mean)
dim(arterial_ec_mean)

#===============================================================================

## Loop for all cell types

library(dplyr)
library(readr)
library(purrr)

pseudo_dir <- "/home/tooba/Skin_Cell_Model/pseudobulk/"

cell_types <- c(
  "Arterial_EC", "Basal_KC", "Bulb", "Bulge",
  "Capillary_EC", "Coil", "Cornified_KC", "Duct",
  "Fibro_A", "Fibro_B", "Fibro_C", "Fibro_D", "Fibro_diseased",
  "Fibro_E", "Granular_KC", "Infundibulum", "Isthmus", "Lymphatic_EC",
  "Prolif__KC", "SG", "SMC", "Spinous_KC", "Venous_EC"
)

compute_mean_expr <- function(ct) {
  f <- file.path(pseudo_dir, paste0("pseudobulk_", ct, ".tsv"))
  df <- read_tsv(f, show_col_types = FALSE)
  donor_cols <- setdiff(colnames(df), "gene_name")
  
  df %>%
    rowwise() %>%
    mutate(!!ct := mean(c_across(all_of(donor_cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    select(gene_name, !!ct)
}

# Compute mean_expression for all 23 cell types
mean_expr_list <- map(cell_types, compute_mean_expr)

# Join into one genes x cell_types matrix
expr_matrix <- reduce(mean_expr_list, full_join, by = "gene_name")

dim(expr_matrix)
head(expr_matrix)


#--------------------------------------
# 1. Confirm all 23 cell type names came through correctly
colnames(expr_matrix)

# 2. Compare against your expected list, to catch typos/mismatches
setdiff(cell_types, colnames(expr_matrix))  # cell types in your list but missing from expr_matrix
setdiff(colnames(expr_matrix), c("gene_name", cell_types))  # anything unexpected

# 3. Confirm each cell type column actually has real, non-degenerate data
# (not all zeros, not all NA)
summary(expr_matrix[ , -1])  # -1 drops gene_name column

# 4. Explicit check: any column that's entirely NA or entirely 0?
sapply(expr_matrix[ , -1], function(x) all(is.na(x)))
sapply(expr_matrix[ , -1], function(x) all(x == 0, na.rm = TRUE))

#---------------------------------------------------------------------

# --------------------------------------------------------------------
# Log-transform mean expression values before correlating cell types.
#
# Why: raw mean_expression has an extreme range within each cell type
# (e.g. Basal_KC ranges from 0 to 112,470). Pearson correlation is
# driven by variance, so a few extremely high-expression genes would
# dominate the result and drown out the pattern across all other
# genes. Log-transforming compresses that range so a 10x change in a
# low-expression gene counts the same as a 10x change in a
# high-expression gene -- which better reflects real biological
# similarity between cell types.
# --------------------------------------------------------------------

# Copy expr_matrix so the original untransformed values are preserved
expr_log <- expr_matrix

# Apply log2(x + 1) to every cell type column (everything except
# gene_name, which is column 1).
# The "+1" avoids log(0), which is undefined -- genes with 0 expression
# become log2(1) = 0 instead of -Inf.
expr_log[ , -1] <- log2(expr_log[ , -1] + 1)

# --------------------------------------------------------------------
# Compute the 23x23 pairwise correlation matrix across cell types.
#
# use = "pairwise.complete.obs": for each pair of cell types, use only
# the genes where BOTH cell types have a non-NA value, rather than
# dropping a gene from the whole matrix if it's missing in just one
# cell type. This handles the NAs introduced by full_join earlier
# without needing a separate filtering step.
# --------------------------------------------------------------------
cor_matrix <- cor(expr_log[ , -1], use = "pairwise.complete.obs")

# Round for easier reading when printed
round(cor_matrix, 2)


# --------------------------------------------------------------------
# Improved heatmap: larger canvas so labels aren't clipped, better
# color contrast since all values sit in the 0.6-1.0 range (a full
# 0-1 color scale wastes most of its range on values you never see),
# and a cleaner border.
# --------------------------------------------------------------------

library(pheatmap)
library(RColorBrewer)

png(
  "/home/tooba/Skin_Cell_Model/celltype_correlation_heatmap.png",
  width = 1800, height = 1500, res = 150   # wider canvas so row labels aren't cut off
)

# Custom color breaks: stretch the color scale across 0.6-1.0 instead
# of 0-1, since that's the actual range of your data. This makes
# differences between e.g. 0.85 and 0.95 visually distinguishable
# instead of both looking like "generic red."
breaks <- seq(0.6, 1.0, length.out = 101)
colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

pheatmap(
  cor_matrix,
  color = colors,
  breaks = breaks,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Cell-type correlation (log2 mean expression)",
  fontsize = 11,            # overall base font size (axis labels, title)
  fontsize_number = 6.5,    # slightly bigger than before, still fits
  cellwidth = 24,
  cellheight = 24,
  border_color = "grey60",  # subtle cell borders instead of none
  angle_col = 45            # readable diagonal column labels
)

dev.off()


##===========================================================================

#Now I want find correlation between cell types using percentiles instead of mean expression
#input files here are processed_ibGs

library(dplyr)
library(readr)
library(purrr)
library(pheatmap)
library(viridis)

processed_dir <- "/home/tooba/Skin_Cell_Model/processed_IbG"

cell_types <- c(
  "Arterial_EC", "Basal_KC", "Bulb", "Bulge",
  "Capillary_EC", "Coil", "Cornified_KC", "Duct",
  "Fibro_A", "Fibro_B", "Fibro_C", "Fibro_D", "Fibro_diseased",
  "Fibro_E", "Granular_KC", "Infundibulum", "Isthmus", "Lymphatic_EC",
  "Prolif__KC", "SG", "SMC", "Spinous_KC", "Venous_EC"
)

# Read gene_name + percentile from one cell type's processed_IbG file
get_percentile <- function(ct) {
  f <- file.path(processed_dir, paste0("processed_IbG_pseudobulk_", ct, ".csv"))
  df <- read_csv(f, show_col_types = FALSE)
  df %>%
    select(gene_name, percentile) %>%
    rename(!!ct := percentile)
}

# Build genes x cell_types matrix of percentiles
pct_list <- map(cell_types, get_percentile)
pct_matrix <- reduce(pct_list, full_join, by = "gene_name")

dim(pct_matrix)
head(pct_matrix)

# Correlate directly -- no log transform needed, percentiles are
# already rank-normalized and bounded 0-1
cor_matrix_pct <- cor(pct_matrix[ , -1], use = "pairwise.complete.obs")

# Heatmap
png(
  "/home/tooba/Skin_Cell_Model/celltype_correlation_heatmap_percentile.png",
  width = 1800, height = 1500, res = 150
)

breaks <- seq(min(cor_matrix_pct), 1.0, length.out = 101)
colors <- viridis(100)

pheatmap(
  cor_matrix_pct,
  color = colors,
  breaks = breaks,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Cell-type correlation (percentile)",
  fontsize = 11,
  fontsize_number = 6.5,
  cellwidth = 24,
  cellheight = 24,
  border_color = "grey60",
  angle_col = 45,
  cutree_rows = 4,
  cutree_cols = 4
)

dev.off()



##================================================================

##to view color pallettes

library(RColorBrewer)
display.brewer.all()          # shows every ColorBrewer palette at once

library(viridis)
# quick visual comparison of the viridis family
n <- 20
par(mfrow = c(4, 1), mar = c(1, 1, 2, 1))
image(1:n, 1, as.matrix(1:n), col = viridis(n), main = "viridis")
image(1:n, 1, as.matrix(1:n), col = magma(n), main = "magma")
image(1:n, 1, as.matrix(1:n), col = plasma(n), main = "plasma")
image(1:n, 1, as.matrix(1:n), col = inferno(n), main = "inferno")

