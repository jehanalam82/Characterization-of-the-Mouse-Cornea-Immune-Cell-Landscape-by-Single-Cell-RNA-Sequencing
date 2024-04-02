###https://cole-trapnell-lab.github.io/monocle3/
# library
library(dplyr)
library(Matrix)
library(ggplot2)
library(reticulate)
library(monocle3)
library(SeuratWrappers)
library(ggridges)
# Set input and output folders
input_folder <- "H:/Manuscripts/2024_K_SCRNA-Methods/Combined_Analysis"
output_folder <- "H:/Manuscripts/2024_K_SCRNA-Methods"
# Specify RDS file path
finrds <- file.path(input_folder, "MFM_Minus_Stromal_legend.rds")
# Create output directory
dirOut <- file.path(output_folder, paste0("out_monocle3_", tools::file_path_sans_ext(basename(finrds))))
dir.create(dirOut, showWarnings = FALSE)  # Create silently if it doesn't exist
# Load data
objA <- readRDS(finrds)
head(objA); objA
# Converting seuratobject to celldataset object for Monocle3
cds <- as.cell_data_set(objA)
head(colData(cds))
#cds <- preprocess_cds(cds, num_dim = 50)
#plot_pc_variance_explained(cds)
#cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "ident")
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds,
           genes = selected_genes,
           color_cells_by = "ident",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
###################Order the cells in pseudotime###########
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

################################################################TOP-DEG#########
marker_test_res <- top_markers(cds, group_cells_by="ident", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
rowData(cds)$gene_short_name <- row.names(rowData(cds))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="ident",
                    ordering_type="maximal_on_diag",
                    max.size=3)

#############fOR MORE THAN ONE TOP MARKERS##########
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="ident",
                    ordering_type="cluster_row_col",
                    max.size=3)

######################################################


save_monocle_objects(cds=cds, directory_path='my_cds_objects', comment='K_Male_Female')
