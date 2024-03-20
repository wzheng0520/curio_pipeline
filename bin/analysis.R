#!/usr/bin/env Rscript

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
T <- TRUE

if (length(args)!=3) {
  stop("sample name, path_to_matched_seurato, max number of beads to be analyzed need to be supplied as arguments", call.=FALSE)
} else  {
  sample = args[1]
  path_to_matched_seurato = args[2]
  #path_to_reference = args[1]
  #max_beads = strtoi(args[4])
  max_beads = strtoi(args[3])
}

library(Seurat)

##* FUNCTIONS
# generate a list of variable features by each cluster
gen_markers_bycluster<-function(variable_features_cluster){
  markers_bycluster<-lapply(unique(variable_features_cluster$cluster), function(x){
    genes<-variable_features_cluster[variable_features_cluster$cluster==x,]
    genes.top<-genes[genes$p_val_adj<0.05,]
    genes.top<-genes.top[order(genes.top$avg_log2FC, decreasing = T),]
  });names(markers_bycluster)<-unique(variable_features_cluster$cluster)
  markers_bycluster
}

# organize spatially variable features obtained via moransi
organize_spatial_features_moransi<-function(seurato){
  variable_features_moransi<-seurato@assays$SCT@meta.features
  variable_features_moransi<-variable_features_moransi[,grep("^[Mm]orans",colnames(variable_features_moransi))]
  variable_features_moransi<-variable_features_moransi[!is.na(variable_features_moransi$MoransI_observed),]
  variable_features_moransi<-variable_features_moransi[order(variable_features_moransi$moransi.spatially.variable.rank),]
  variable_features_moransi
}

# add spatial info to seurat object
add_spatial<-function(seurato, matched_barcode){
  beads_ordered <- Cells(seurato)
  seurato[["SPATIAL"]] <- CreateDimReducObject(embeddings = as.matrix(matched_barcode)[beads_ordered,c(1,2)],
                                               key = "SPATIAL_", assay = DefaultAssay(seurato))
  coords = data.frame(x=-matched_barcode[beads_ordered, 2], y=matched_barcode[beads_ordered, 1],
                      row.names=beads_ordered, stringsAsFactors=FALSE)
  seurato@images$slice1 = new(Class = "SlideSeq", assay = "Spatial", key = "slice1_", coordinates = coords)
  seurato
}

# rotate tile to match physical orientation
rotate_whitelist <- function(matched_barcode){
  adjusted_matched_barcode <- data.frame("SPATIAL_1"=-matched_barcode[,2],"SPATIAL_2"=-matched_barcode[,1])
  rownames(adjusted_matched_barcode) <- rownames(matched_barcode)
  adjusted_matched_barcode
}

# quantify genes of interest per base
count_house_keeping_genes <- function(genes_of_interest) {
  tryCatch(
          {
            pct_genes<-colSums(matched_seurato[["RNA"]]@counts[rownames(matched_seurato[["RNA"]]@counts)%in%genes_of_interest,])/colSums(matched_seurato[["RNA"]]@counts)
            return(pct_genes)
          },
           error = function(c) {
             msg=conditionMessage(c)
             print(paste0(msg))
             pct_genes<-rep(0, length(Cells(matched_seurato)))
             names(pct_genes)<-Cells(matched_seurato)
             return(pct_genes)
           },
           warning = function(c) "warning",
           message = function(c) "message"
  )
}

# import genes of interest
import_house_keeping_genes <- function(path_to_genes) {
  tryCatch(
          {
            genes<-read.table(path_to_genes, header=F)[,1]
            return(genes)
          },
           error = function(c) {
             msg=conditionMessage(c)
             print(paste0(msg))
             genes=c()
             return(genes)
           },
           warning = function(c) "warning",
           message = function(c) "message"
  )
}

##############################################################
# DATA IMPORT
##############################################################
##* DATA IMPORT
message("Reading seurat rds")
matched_seurato <- readRDS(path_to_matched_seurato)

num_beads <- length(Cells(matched_seurato))

if(num_beads > max_beads){
   message("More than ", max_beads, " beads were recovered. Subsampling to ",max_beads, " beads for analysis.")
   set.seed(123)
   matched_seurato <- matched_seurato[, sample(colnames(matched_seurato), size = max_beads, replace=F)]
}

print("Matched seurat")
print(head(matched_seurato))

##############################################################
# PRIMARY
##############################################################
message("Adding log10_nCount_RNA")
matched_seurato <- AddMetaData(matched_seurato, log(matched_seurato$nCount_RNA,base=10), col.name = c("log10_nCount_RNA"))
message("Adding log10_nFeature_RNA")
matched_seurato <- AddMetaData(matched_seurato, log(matched_seurato$nFeature_RNA,base=10), col.name = c("log10_nFeature_RNA"))
message("Running SCTransform")
matched_seurato <- SCTransform(
                        matched_seurato,
                        assay = "RNA",
#                         ncells = 1000,
#                         vars.to.regress = "percent.mt",
                        verbose = TRUE,
                        conserve.memory = TRUE
                    )
message("Running PCA")
matched_seurato <- RunPCA(matched_seurato)

##############################################################
# CLUSTERING
##############################################################
message("Clustering")
n_pc <- ncol(Loadings(matched_seurato, reduction = "pca"))
n_dim <- min(n_pc, 30)
matched_seurato$log_umi <- matched_seurato$log10_nCount_RNA
matched_seurato <- RunUMAP(matched_seurato, dims = 1:n_dim)
matched_seurato <- FindNeighbors(matched_seurato, dims = 1:n_dim)
matched_seurato <- FindClusters(matched_seurato, resolution = 0.2)
variable_features_cluster <- FindAllMarkers(matched_seurato, assay = "SCT", only.pos =T)
message("Complete FinalAllMarkers")

## order cluster ID numerically
DefaultAssay(matched_seurato)<-"SCT"
matched_seurato$SCT_snn_res.0.2<-factor(matched_seurato$SCT_snn_res.0.2, levels=(sort(as.numeric(levels(matched_seurato$SCT_snn_res.0.2)))))

## identify spatially variable genes and cleanup
library(Rfast2)
matched_seurato <- FindSpatiallyVariableFeatures(matched_seurato, assay = "SCT", slot = "scale.data", features = VariableFeatures(matched_seurato)[1:200], selection.method = "moransi", x.cuts = 100, y.cuts = 100, verbose = TRUE, nfeatures=200)
variable_features_moransi<-organize_spatial_features_moransi(matched_seurato)

##* OUTPUT
saveRDS(matched_seurato, paste0(sample, "_","seurat.rds"))
sceasy::convertFormat(matched_seurato, from="seurat", to="anndata", outFile=paste0(sample, "_", "anndata.h5ad"))
write.table(matched_seurato$seurat_clusters, paste0(sample, "_","cluster_assignment.txt"), quote = F, sep="\t", col.names = F)
write.table(variable_features_cluster, paste0(sample, "_", "variable_features_clusters.txt"), quote = F, sep="\t")
write.table(variable_features_moransi, paste0(sample, "_","variable_features_spatial_moransi.txt"), quote = F, sep="\t")
