#heatmap of TSS
#libraries
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(peakSeason)

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"
input.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/6_tracks_normalized"
#load data
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#Set specifications 
cutoff <- 0.01 #set FDR cutoff
l2fc <- 2##set logfold2 cutoff

#annotations
#anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
#anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
repeats <- fread("/omics/groups/OE0219/internal/repguide/data/repeats/hg19_repeats.txt")
responsive_TSS <- readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/responsive_TSSs.RDS")
