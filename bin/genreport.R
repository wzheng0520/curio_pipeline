#!/usr/bin/env Rscript

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
print(args)

if (length(args)!=14) {
  stop("14 arguments need to be supplied as arguments", call.=FALSE)
} else  {
  sample = args[1]	
  file_whitelist <- args[2]
  beadbarcode_file <- strsplit(file_whitelist, split="/")[[1]][length(strsplit(file_whitelist, split="/")[[1]])]
  beadbarcode <- gsub("_BeadBarcodes.txt", "", beadbarcode_file)

  rmarkdown::render( 
	input  = args[3], 
	output_file = paste0(sample, "_", "Report.html"),
       	output_dir = getwd(),
	intermediates_dir = getwd(),
	knit_root_dir = getwd(),	
	params = list( 
	    sample = sample,
	    beadbarcode = beadbarcode,
	    file_metrics = args[4],
	    file_matched_seurato = args[5],
	    file_variable_features_clusters = args[6],
	    file_variable_features_spatial_moransi = args[7],
	    file_read1_properStructure = args[8],
	    file_read1_improperStructure = args[9],
	    file_numReads_perBD = args[10],
	    file_numReads_perBB = args[11],
	    file_BB_BM_BD = args[12],
	    reference = args[13],
	    seekerexplorer = args[14]
	)
 )
}
