#!/usr/bin/env Rscript

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
print(length(args))

if (args[1] == "Seurat_default") {
  if (length(args) == 4) {
    algorithm = args[1]
    sample = args[2]
    path_to_seurat_object = args[3]
    path_to_spatial_var_score = args[4]
  } else { 
    stop("algorithm, sample_name, path to seurat_object, path to spatial_var_scores need to be supplied as arguments",call.=FALSE)
  }
} else if (args[1] == "RCTD") {
    if (length(args) == 5) {
      algorithm = args[1]
      sample = args[2]
      path_to_seurat_object = args[3]
      path_to_spatial_var_score = args[4]
      path_to_reference = args[5]
   } else {
     stop("algorithm, sample_name, path to seurat_object, path to spatial_var_scores, path to ref(RCTD)/path to seurat_cluster_assignments(Seurat_default) need to be supplied as arguments",call.=FALSE)
   }
}
# write the above only for seurat
library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(utils)
library(arrow)
library(stats)
library(dplyr)

##* FUNCTIONS

# generate output files for Seurat_default or RCTD
process_Seeker_data <- function(algorithm, sample_name, seurat_object, spatial_var_scores, reference = NULL, UMI_min = 0, UMI_min_sigma = 300, gene_count_min = 0, max_cores = 8) {
  if (algorithm == "Seurat_default") {
    cluster_assignments <- data.frame(sobj[['seurat_clusters']])
    df_Seurat_parquet <- create_Parquet_df(count_matrix, coords, sobj, moransi_scores, UMI_min, gene_count_min, cluster_assignments)
    return(list(df_Seurat_parquet))
  } else if (algorithm == "RCTD") {
    spatial_RNA <- SpatialRNA(coords, count_matrix, nUMI)
    RCTD <- create.RCTD(spatial_RNA, readRDS(path_to_reference), max_cores, CELL_MIN_INSTANCE = 0, UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    RCTD_clusters <- as.data.frame(RCTD@results$results_df["first_type"])
    confidence_levels <- as.data.frame(RCTD@results$results_df["spot_class"])
    df_RCTD_parquet <- create_Parquet_df(count_matrix, coords, sobj, moransi_scores,UMI_min, gene_count_min, RCTD_clusters, confidence_levels)
    df_annotation <- create_annotate_file(as.data.frame(RCTD@results$results_df))
    df_seurat <- data.frame(RCTD@results$results_df["first_type"], RCTD@results$results_df["second_type"], RCTD@results$results_df["spot_class"])
    sobj  <- AddMetaData(object = sobj , metadata = df_seurat[Cells(sobj),], col.name=colnames(df_seurat))
    plot_spatial<- DimPlot(sobj, reduction = 'SPATIAL', group.by = 'first_type')&NoAxes()&theme(plot.title = element_blank())&coord_fixed()
    return(list(RCTD, df_RCTD_parquet, df_annotation, sobj, plot_spatial))
  }
  print("Seeker data has been processed. Find the results in the output directory.")
}

# create parquet file
create_Parquet_df <- function(counts, spatial_coordinates, seurat_object, spatial_var_scores, UMI_min, gene_count_min, cluster_assignments, cluster_assignment_confidence_levels = NULL) {
    count_matrix <- as.data.frame(t(as.matrix(counts)))
    UMI_counts <- rowSums(count_matrix)
    gene_counts <- rowSums(count_matrix != 0)
    count_matrix["UMI_count"] <- UMI_counts
    count_matrix["gene_count"] <- gene_counts
    colnames(spatial_coordinates) <- c("SPATIAL_1", "SPATIAL_2")
    UMAP_coordinates <- as.data.frame(seurat_object[["umap"]]@cell.embeddings)
    colnames(cluster_assignments) <- c("cluster")
    if (is.null(cluster_assignment_confidence_levels)) {
      df <- merge(cluster_assignments, cbind(spatial_coordinates,
      UMAP_coordinates, count_matrix), by = 0)
      colnames(df) <- c("barcode", colnames(df)[2:length(df)])
      df <- subset(df, UMI_count >= UMI_min & gene_count >= gene_count_min)
      total_abundance <- colSums(df[7:length(df)])
      cv <- sapply(df[7:length(df)], stats::sd) / colMeans(df[7:length(df)])
    } else {
      colnames(cluster_assignment_confidence_levels) <- c("confidence")
      df <- merge(cbind(cluster_assignments, cluster_assignment_confidence_levels), cbind(spatial_coordinates, UMAP_coordinates, count_matrix), by = 0)
      colnames(df) <- c("barcode", colnames(df)[2:length(df)])
      df <- subset(df, UMI_count >= UMI_min & gene_count >= gene_count_min)
      total_abundance <- colSums(df[8:length(df)])
      cv <- sapply(df[8:length(df)], stats::sd) / colMeans(df[8:length(df)])
    }
    moransi_scores <- unlist(spatial_var_scores)
    names(moransi_scores) <- rownames(spatial_var_scores)
    df <- dplyr::bind_rows(df, total_abundance, cv, moransi_scores)
    return(df[order(df$cluster), ])
}

# create annotation file
create_annotate_file <- function(df) {
  df_anno <- data.frame(rownames(df), df["first_type"], df["second_type"], df["spot_class"])
  colnames(df_anno) <- c("Barcode", "First_Cell_Type", "Second_Cell_Type", "Confidence_Level")
  #write.table(df_anno, file = file_name, sep="\t", row.names = FALSE,
              #col.names = TRUE, quote = FALSE)
  return(df_anno)
}

##* DATA IMPORT
#cluster_assignments <- utils::read.table(path_to_reference, row.names = 1)
sobj <- readRDS(path_to_seurat_object) #load seurat object
total_beads <- ncol(sobj)
if (total_beads > 100000) {
   set.seed(1) 
   sobj <- sobj[, sample(colnames(sobj), size = total_beads/10, replace=F)]
}
count_matrix <- (sobj[["RNA"]]@counts) #get matrix using seurat object
nUMI <- colSums(count_matrix) #get nUMI - list of total counts
coords <- as.data.frame(sobj[["SPATIAL"]]@cell.embeddings) #get spatial_coordinates using seurat object
moransi_scores <- utils::read.table(path_to_spatial_var_score, row.names = 1)[, 1, drop = FALSE] #read spatial_var_scores file


## ANALYSIS
process_data <- process_Seeker_data (algorithm = algorithm, sample_name = sample, seurat_object = path_to_seurat_object, spatial_var_scores = path_to_spatial_var_score,
                     reference = path_to_reference, gene_count_min = 0, max_cores = 8)

##* OUTPUT
if (algorithm == "Seurat_default") {
   arrow::write_parquet(process_data[[1]], paste0(sample, "_", "Seurat_default", "_", "parquet.gzip"))
} else if (algorithm == "RCTD") { 
   arrow::write_parquet(process_data[[2]], paste0(sample, "_", "RCTD", "_", "parquet.gzip")) 
   write.table(process_data[[3]], paste0(sample, "_", "RCTD", "_", "annotation.txt"), sep="\t", row.names = FALSE,col.names = TRUE, quote = FALSE)
   saveRDS(process_data[[1]], paste0(sample, "_", "RCTD", ".rds" ))
   saveRDS(process_data[[4]], paste0(sample, "_", "RCTD", "_", "seurat.rds"))
   png(paste0(sample, "_", "RCTD", "_", "spatial.png"), units="in", width=12, height=20, res=1200)
   plot(process_data[[5]])
   dev.off()
}
















