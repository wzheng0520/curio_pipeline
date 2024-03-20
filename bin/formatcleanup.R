#!/usr/bin/env Rscript

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("sample name, expression matrix, numReads_perBM, path to whitelist, path to ref, minUMI_perBead, bead_chunk_size need to be supplied as arguments", call.=FALSE)
} else  {
  sample = args[1]
  expression_matrix = args[2]
  numReads_perBM = args[3]
  path_to_whitelist = args[4]
  path_to_reference = args[5]
  minUMI_perBead = strtoi(args[6])
  bead_chunk_size = strtoi(args[7])
}

library(Seurat)
library(Matrix)
library(data.table)
library(arrow)
library(dplyr)

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

##* DATA IMPORT
# split columns ids in to chunks
DS_numReads_perBM <- arrow::open_dataset(sources = numReads_perBM)
num_col <- nrow(DS_numReads_perBM)+1
print(paste0("Number of columns: ", num_col))

col_splits <- split(2:num_col,ceiling(seq_along(2:num_col) / bead_chunk_size))

message("read in parquet file")
DS <- arrow::open_dataset(sources = expression_matrix)

# read in genes
genes <- DS %>% select("gene")
genes <- Scanner$create(genes)
genes <- genes$ToTable()
genes <- as.character(genes$gene)
print(paste0("Number of genes: ", length(genes)))

message("----------------------------------------------------------")
message("save table to seurat objects in chunks")
list_seurat_matched_dt <- lapply(col_splits, function(x){
  DS_t <- DS %>% select(all_of(x))
  SO <- Scanner$create(DS_t)
  AT <- SO$ToTable()
  matched_dt <- as.data.frame(AT)
  print(dim(matched_dt))
  barcodes_toretain <- which(colSums(matched_dt) >= minUMI_perBead)
  matched_dt <- matched_dt[,barcodes_toretain,drop=FALSE]
  matched_dt <- data.frame(matched_dt)
  print(dim(matched_dt))

  if(dim(matched_dt)[2]>0){
  	rownames(matched_dt) <- genes
 	matched_seurat_dt <- CreateSeuratObject(counts = matched_dt, project = sample)
  	matched_seurat_dt
  }else{
	NULL
  }
})

list_seurat_matched_dt<-list_seurat_matched_dt[!sapply(list_seurat_matched_dt,is.null)]
message("complete save table to seurat objects in chunks")
message("----------------------------------------------------------")

message("----------------------------------------------------------")
message("merge seurat objects")
if(length(list_seurat_matched_dt)==1){
  	matched_seurato <- list_seurat_matched_dt[[1]]
}else{
  	#matched_seurato <- merge(list_seurat_matched_dt[[1]], y = list_seurat_matched_dt[-1], add.cell.ids = rep("", length(list_seurat_matched_dt)), project = sample)
  	matched_seurato <- merge(list_seurat_matched_dt[[1]], y = list_seurat_matched_dt[-1], add.cell.ids = NULL, project = sample)
	print(length(Cells(matched_seurato)))
	#matched_seurato <- RenameCells(matched_seurato, new.names=sub('^_', '', Cells(matched_seurato)))  #remove underscores from start of cell names
}
rm("list_seurat_matched_dt")
message("complete merge seurat objects")
message("----------------------------------------------------------")

whitelist_coord <- read.table(path_to_whitelist, header=FALSE, row.names = 1)
matched_barcode <- whitelist_coord[rownames(whitelist_coord)%in%Cells(matched_seurato),]
colnames(matched_barcode) <- c("SPATIAL_1","SPATIAL_2")
rm("whitelist_coord")

matched_barcode_adjusted <- rotate_whitelist(matched_barcode)

##* ADD METADATA
# spatial
matched_seurato <- add_spatial(matched_seurato, matched_barcode_adjusted)
rm("matched_barcode_adjusted")

### nReads perCell
SO <- Scanner$create(DS_numReads_perBM)
AT <- SO$ToTable()

matched.numReads_perCell <- data.frame(AT)
rm(DS_numReads_perBM,SO,AT)

print(head(matched.numReads_perCell))
row.names(matched.numReads_perCell)<-matched.numReads_perCell$BM
# matched.numReads_perCell<-matched.numReads_perCell[,c(1), drop = FALSE]
matched.numReads_perCell$BM_count = matched.numReads_perCell$BM_count


print(head(matched.numReads_perCell))
matched.numReads_perCell$log <- log(matched.numReads_perCell[,"BM_count"], base=10)
matched_seurato <- AddMetaData(matched_seurato, matched.numReads_perCell[rownames(matched.numReads_perCell)%in%names(matched_seurato$orig.ident), c("BM_count","log")], col.name = c("numReads", "log10_numReads"))
rm("matched.numReads_perCell")

# log transform nCount and nFeature
matched_seurato <- AddMetaData(matched_seurato, log(matched_seurato$nCount_RNA,base=10), col.name = c("log10_nCount_RNA"))
matched_seurato <- AddMetaData(matched_seurato, log(matched_seurato$nFeature_RNA,base=10), col.name = c("log10_nFeature_RNA"))

# ribosomal, mitochondial genes
genes_rrna <- import_house_keeping_genes(file.path(path_to_reference, "rRNA_genes.txt"))
genes_mito <- import_house_keeping_genes(file.path(path_to_reference, "mt_genes.txt"))
genes_mito <- genes_mito[!genes_mito %in% intersect(genes_rrna, genes_mito)]

# For handling concatenated genes and transcripts
genes_rrna <- sub("|", "-", genes_rrna, fixed=TRUE)
genes_rrna <- sub("_", "-", genes_rrna, fixed=TRUE)
genes_mito <- sub("|", "-", genes_mito, fixed=TRUE)
genes_mito <- sub("_", "-", genes_mito, fixed=TRUE)

# species: "Mus_musculus",""Homo_sapiens","Gallus_gallus","Danio_rerio"
if(grepl("Mus_musculus", path_to_reference)){
        genes_rpro <- rownames(matched_seurato[["RNA"]]@counts)[grep("^Rp[sl]", rownames(matched_seurato[["RNA"]]@counts))]
}else if(grepl("Homo_sapiens", path_to_reference)){
        genes_rpro <- rownames(matched_seurato[["RNA"]]@counts)[grep("^RP[SL]", rownames(matched_seurato[["RNA"]]@counts))]
}else if(grepl("Gallus_gallus", path_to_reference)){
        genes_rpro <- rownames(matched_seurato[["RNA"]]@counts)[grep("^RP[SL]", rownames(matched_seurato[["RNA"]]@counts))]
}else if(grepl("Zea_mays|Danio_rerio", path_to_reference)){
	genes_rpro <- rownames(matched_seurato[["RNA"]]@counts)[grep("^rp[sl]", rownames(matched_seurato[["RNA"]]@counts))]
}else{
	genes_rpro <- c()
        message("Reference species unrecognized, check reference to make sure MT and RRNA genes are defined ...")
}

pct_rrna <- count_house_keeping_genes(genes_rrna)
pct_rpro <- count_house_keeping_genes(genes_rpro)
pct_mito <- count_house_keeping_genes(genes_mito)

matched_seurato <- AddMetaData(matched_seurato, pct_rrna, col.name = c("pct_rrna"))
matched_seurato <- AddMetaData(matched_seurato, pct_rpro, col.name = c("percent.rp"))
matched_seurato <- AddMetaData(matched_seurato, pct_mito, col.name = c("percent.mt"))

rm(list=c("genes_rrna","genes_mito","genes_rpro","pct_rrna", "pct_rpro", "pct_mito"))
message("----------------------------------------------------------")
message(head(matched_seurato))
message("----------------------------------------------------------")
message("Saving seurat object")

# output
saveRDS(matched_seurato, paste0(sample, "_","seurat.rds"))
writeMM(obj=matched_seurato[["RNA"]]@counts, paste0(sample, "_", "MoleculesPerMatchedBead.mtx"))
write(x = rownames(matched_seurato[["RNA"]]@counts), file = paste0(sample, "_", "genes.tsv"))
write(x = colnames(matched_seurato[["RNA"]]@counts), file = paste0(sample, "_", "barcodes.tsv"))
write.csv(matched_barcode, paste0(sample, "_", "MatchedBeadLocation.csv"), quote = F)

