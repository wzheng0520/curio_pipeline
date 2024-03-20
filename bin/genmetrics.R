#!/usr/bin/env Rscript

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)!=9) {
  stop("sample ID, alignmentqc_proper, whitelist, matched_seurat_object, numReads_perBD, numReads_perXS, featureCount report, report format, scripts folder need to be supplied as an argument", call.=FALSE)
} else {
  sample = args[1]
  file_whitelist = args[2]
  file_readsummary = args[3]
  file_matched_seurato = args[4]
  file_numReads_perXS = args[5] #matched
  file_FC_log =  args[6]
  file_numReads_perBD = args[7] #matched
  report = args[8]
  scripts_folder = args[9]
}

library(Seurat)

# read summaries
summary_r1 <- read.table(file_readsummary, header=F, sep="\t", stringsAsFactors = F, row.names = 1)
matched_seurato <- readRDS(file_matched_seurato)
numReads_perXS <- read.table(file_numReads_perXS, header=F, sep="\t", stringsAsFactors = F, row.names = 1)
FC_log <- read.table(file_FC_log, header=T, sep="\t", stringsAsFactors = F, row.names = 1)
whitelist <- read.table(file_whitelist, header=FALSE, row.names = 1)

numBarcodes_on_whitelist <- nrow(whitelist)
rm("whitelist")
numBarcodes_matched <- length(Cells(matched_seurato))

message("read quality")
df_read_quality<-data.frame(
  "Total_raw_fastq_readpairs"=c(summary_r1["total", 1]),
  "Readpairs_proper_structure"=c(summary_r1["proper_structure", 1]),
  "Readpairs_with_adapter"=c(summary_r1["total", 1] - summary_r1["keep_illumina", 1]),
  "Percent_readpairs_proper_structure"=c(100*round(summary_r1["proper_structure", 1]/summary_r1["total", 1], 4)),
  "Percent_readpairs_with_adapter"=c(100*round((1 - summary_r1["keep_illumina", 1]/summary_r1["total", 1]), 4)),  
  "Percent_adapter_free_readpairs_proper_structure"=c(100*round(summary_r1["proper_structure", 1]/summary_r1["keep_illumina", 1], 4)),
  "Total_reads_in_genes"=c(FC_log["Assigned", 1]), #adjusted
  "Percent_reads_in_genes"=c(100*round(FC_log["Assigned", 1]/summary_r1["total", 1], 4)) #adjusted
);row.names(df_read_quality)<-"Values"

message("bead matching quality")
df_bead_matching_quality<-data.frame(
  "Total_barcodes_on_substrate"=c(numBarcodes_on_whitelist),
  "Total_barcodes_recovered"=c(numBarcodes_matched),
  "Percent_barcodes_recovered"=c(100*round(numBarcodes_matched/numBarcodes_on_whitelist, 4)),
  "Reads_matched_to_barcodes"=c(summary_r1["r1_proper_structure_matched", 1]),
  "Percent_proper_reads_matched_to_barcodes"=c(100*round(summary_r1["r1_proper_structure_matched", 1]/summary_r1["proper_structure", 1],4))
);row.names(df_bead_matching_quality)<-"Values"

message("gene alignment quality") 
df_gene_alignment_quality<-data.frame(				      
  "Reads_matched_to_barcodes_and_genic" = c(numReads_perXS["Assigned", 1]),
  "Percent_proper_reads_matched_to_barcodes_and_genic"=c(100*round(numReads_perXS["Assigned", 1]/summary_r1["r1_proper_structure_matched", 1],4))
);row.names(df_gene_alignment_quality)<-"Values"

message("overall sample quality ")
df_overall_sample_quality<-data.frame(				       
  "Percent_raw_reads_useful_matched_and_genic"=c(100*round(numReads_perXS["Assigned", 1]/summary_r1["total", 1], 4)),
  "Median_reads_perUMI"=c(round(median(matched_seurato$numReads)/median(matched_seurato$nCount_RNA),2)),
  "Total_genes_in_matched_bead_barcodes"=c(sum(rowSums(matched_seurato[["RNA"]]@counts)>0)),
  "Total_UMI_in_matched_bead_barcodes"=c(sum(matched_seurato$nCount_RNA)),
	
  "Mean_reads_perBead"=c(round(mean(matched_seurato$numReads),2)),
  "Mean_UMI_perBead"=c(round(mean(matched_seurato$nCount_RNA),2)),
  "Mean_gene_perBead"=c(round(mean(matched_seurato$nFeature_RNA),2)),
  "Mean_reads_perUMI"=c(round(mean(matched_seurato$numReads)/mean(matched_seurato$nCount_RNA),2)), 
  
  "Median_reads_perBead"=c(round(median(matched_seurato$numReads),2)),
  "Median_UMI_perBead"=c(round(median(matched_seurato$nCount_RNA),2)),
  "Median_gene_perBead"=c(round(median(matched_seurato$nFeature_RNA),2)),
  
  "Top10_percent_reads_perBead"=c(round(quantile(matched_seurato$numReads, 0.9),2)),
  "Top10_percent_UMI_perBead"=c(round(quantile(matched_seurato$nCount_RNA, 0.9),2)),
  "Top10_percent_gene_perBead"=c(round(quantile(matched_seurato$nFeature_RNA, 0.9),2)),

  "Median_percent_mitochondrial_UMI_perBead"=c(round(median(matched_seurato$percent.mt),4)*100),
  "Median_percent_ribosomalProtein_UMI_perBead"=c(round(median(matched_seurato$percent.rp),4)*100),
  "Median_percent_ribosomalRNA_UMI_perBead"=c(round(median(matched_seurato$pct_rrna),4)*100)
);row.names(df_overall_sample_quality)<-"Values"

message("Compiling metrics and output ...")
list_quality<-list("read_quality"=df_read_quality,
                   "bead_matching_quality"=df_bead_matching_quality,
                   "gene_alignment_quality"=df_gene_alignment_quality,
                   "overall_sample_quality"=df_overall_sample_quality)

list_quality<-lapply(1:length(list_quality), function(x){
  list_quality[[x]]<-as.data.frame(t(list_quality[[x]]))
  list_quality[[x]]$Category<-names(list_quality)[x]
  list_quality[[x]]
})

matrix_quality<-do.call("rbind", list_quality)

if(report == "d"){
	numReads_perBD = read.table(file_numReads_perBD, header=F, sep="\t", stringsAsFactors = F)
	numReads_perfectly_matched = numReads_perBD[numReads_perBD[,1]==0,2]
	rm("numReads_perBD")
	matrix_quality = rbind(matrix_quality, c(numReads_perfectly_matched, "bead_matching_quality"))
	rownames(matrix_quality)[nrow(matrix_quality)] = "Reads_perfectly_matched_to_barcodes"
	print(matrix_quality)
}

if (report == "d") {
  	source(file.path(scripts_folder,"peak_detector",'py_caller.R'))
	peak_count <- as.numeric(call_py_script(paste0(scripts_folder, '/peak_detector/py_script/r_transition.py'), matched_seurato$log10_nCount_RNA, paste0('./temp_', sample, '.csv')))
	print(paste0("Number of peak detected: ", peak_count)) 
	if (peak_count>1) {
        	d<-density(matched_seurato$log10_nCount_RNA)
        	cutoff_log10_UMI<-optimize(approxfun(d$x,d$y),interval=c(1,quantile(matched_seurato$log10_nCount_RNA, 0.9)))$minimum
        	names_high<-names(matched_seurato$log10_nCount_RNA)[matched_seurato$log10_nCount_RNA>=cutoff_log10_UMI]
        	names_low<-names(matched_seurato$log10_nCount_RNA)[matched_seurato$log10_nCount_RNA<cutoff_log10_UMI]
        	knee_bead <- names(sort(matched_seurato$log10_nCount_RNA[names_high], decreasing = F))[1]

        	df_update<-as.data.frame(t(data.frame("Number_of_peaks"=peak_count, "Total_barcodes_included_by_knee"=length(names_high), "cutoff_log10_UMI"=cutoff_log10_UMI)))
        	df_update$Category<-"overall_library_quality"
        	colnames(df_update)[1]<-"Values"
        	#print(df_update)
        	matrix_quality["Median_reads_perUMI", "Values"]=c(round(median(matched_seurato$numReads[names_high])/median(matched_seurato$nCount_RNA[names_high]),2))
        	matrix_quality["Median_reads_perBead", "Values"]=c(round(median(matched_seurato$numReads[names_high]),2))
        	matrix_quality["Median_UMI_perBead", "Values"]=c(round(median(matched_seurato$nCount_RNA[names_high]),2))
        	matrix_quality["Median_gene_perBead", "Values"]=c(round(median(matched_seurato$nFeature_RNA[names_high]),2))

        	matrix_quality["Top10_percent_reads_perBead", "Values"]=c(round(quantile(matched_seurato$numReads[names_high], 0.9),2))
        	matrix_quality["Top10_percent_UMI_perBead", "Values"]=c(round(quantile(matched_seurato$nCount_RNA[names_high], 0.9),2))
        	matrix_quality["Top10_percent_gene_perBead", "Values"]=c(round(quantile(matched_seurato$nFeature_RNA[names_high], 0.9),2))
	} else {
        	df_update<-as.data.frame(t(data.frame("Number_of_peaks"=peak_count, "Total_barcodes_included_by_knee"=length(matched_seurato$log10_nCount_RNA), "cutoff_log10_UMI"=min(matched_seurato$log10_nCount_RNA))))
        	df_update$Category<-"overall_library_quality"
        	colnames(df_update)[1]<-"Values"
	}
	matrix_quality<-rbind(matrix_quality, df_update)
}

write.csv(matrix_quality, paste0(sample, "_", "Metrics.csv"), quote=F)










