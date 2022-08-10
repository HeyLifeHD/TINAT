#epi heatmaps
#libraries
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(GetoptLong)
library(circlize)
library(RColorBrewer)
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
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
#Set specifications 
cutoff <- 0.01 #set FDR cutoff
l2fc <- 2##set logfold2 cutoff

#parameters
#region
ext <- 1500 
fontSize <- 10
#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
#get regions of interest 
anno_classi$transcript_id <- anno_classi$qry_id
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
DEG_results_list_anno_sub <- lapply(DEG_results_list, function(x){
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    x
})
DEG_results_list_anno_sub_sub <- DEG_results_list_anno_sub$DACandSB939_vs_DMSO
anno_transcript <- anno[anno$type=="transcript",]
tss <- resize(anno_transcript[anno_transcript$transcript_id %in% DEG_results_list_anno_sub_sub$transcript_id, ],1)
tss$class_code_simple <- DEG_results_list_anno_sub_sub$class_code_simple

#rnaseq -> merged tracks
DMSO_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DMSO.normalized.bigWig")
SB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/SB939.normalized.bigWig")
DAC_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DAC.normalized.bigWig")
DACSB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DACSB939.normalized.bigWig")

#cage
DMSO_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DMSO.bigWig")
SB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/SB939.bigWig")
DAC_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DAC.bigWig")
DACSB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DACSB.bigWig")

#methylation
DMSO_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
SB_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150389_H2_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
DAC_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
DACSB_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150387_H2_DAC_plus_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")

#H3K9me3
DMSO_9me3 <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150368_DMSO_H3K9me3_ATCACG_L002_R1_complete_filtered.fastq.gz.bigWig")
SB_9me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150384_SB939_H3K9me3_TTAGGC_L002_R1_complete_filtered.fastq.gz.bigWig")
DAC_9me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150352_DAC_H3K9me3_CGATGT_L002_R1_complete_filtered.fastq.gz.bigWig")
DACSB_9me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150360_DACSB_H3K9me3_ACTTGA_L002_R1_complete_filtered.fastq.gz.bigWig")

#H3K4me3
DMSO_4me3 <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150366_DMSO_H3K4me3_ATCACG_L001_R1_complete_filtered.fastq.gz.bigWig")
SB_4me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150382_SB939_H3K4me3_ACTTGA_L001_R1_complete_filtered.fastq.gz.bigWig")
DAC_4me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150350_DAC_H3K4me3_GCCAAT_L001_R1_complete_filtered.fastq.gz.bigWig")
DACSB_4me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150358_DACSB_H3K4me3_CTTGTA_L001_R1_complete_filtered.fastq.gz.bigWig")

#H3K9ac
DMSO_9ac <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150367_DMSO_H3K9ac_CGATGT_L005_R1_complete_filtered.fastq.gz.bigWig")
SB_9ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150383_SB939_H3K9ac_TGACCA_L005_R1_complete_filtered.fastq.gz.bigWig")
DAC_9ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150351_DAC_H3K9ac_TTAGGC_L005_R1_complete_filtered.fastq.gz.bigWig")
DACSB_9ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150359_DACSB_H3K9ac_GATCAG_L005_R1_complete_filtered.fastq.gz.bigWig")

#H3K27ac
DMSO_27ac <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150362_DMSO_H3K27ac_ATCACG_L008_R1_complete_filtered.fastq.gz.bigWig")
SB_27ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150378_SB939_H3K27ac_ACTTGA_L008_R1_complete_filtered.fastq.gz.bigWig")
DAC_27ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150346_DAC_H3K27ac_GCCAAT_L008_R1_complete_filtered.fastq.gz.bigWig")
DACSB_27ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150354_DACSB_H3K27ac_CTTGTA_L008_R1_complete_filtered.fastq.gz.bigWig")

#H3K14ac
chip.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/"
DMSO_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597665_DMSO-7.normalized", ".bigWig")))
SB_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597647_SB-7.normalized", ".bigWig")))
DAC_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597638_DAC-7.normalized", ".bigWig")))
DACSB_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597656_DAC+SB-7.normalized", ".bigWig")))

#H2BK5ac
DMSO_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597670_DMSO-15.normalized", ".bigWig")))
SB_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597652_SB-15.normalized", ".bigWig")))
DAC_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597643_DAC-15.normalized", ".bigWig")))
DACSB_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597661_DAC+SB-15.normalized", ".bigWig")))

#create bigwig list
bw_chip_list <- list(DMSO_rnaseq=DMSO_rnaseq, SB_rnaseq=SB_rnaseq, DAC_rnaseq=DAC_rnaseq, DACSB_rnaseq=DACSB_rnaseq,
    DMSO_cage=DMSO_cage, SB_cage=SB_cage, DAC_cage=DAC_cage, DACSB_cage=DACSB_cage,
    DMSO_9me3=DMSO_9me3, SB_9me3=SB_9me3, DAC_9me3=DAC_9me3, DACSB_9me3=DACSB_9me3,
    DMSO_4me3=DMSO_4me3,SB_4me3=SB_4me3, DAC_4me3=DAC_4me3, DACSB_4me3=DACSB_4me3,
    DMSO_9ac=DMSO_9ac, SB_9ac=SB_9ac, DAC_9ac=DAC_9ac, DACSB_9ac=DACSB_9ac, 
    DMSO_27ac=DMSO_27ac, SB_27ac=SB_27ac, DAC_27ac=DAC_27ac, DACSB_27ac=DACSB_27ac,
    DMSO_14ac=DMSO_14ac, SB_14ac=SB_14ac, DAC_14ac=DAC_14ac, DACSB_14ac=DACSB_14ac,
    DMSO_5ac=DMSO_5ac, SB_5ac=SB_5ac, DAC_5ac=DAC_5ac, DACSB_5ac=DACSB_5ac
     )
bw_meth_list <- list(DMSO_meth=DMSO_meth, SB_meth=SB_meth, DAC_meth=DAC_meth, DACSB_meth=DACSB_meth) 
bw_meth_list <- lapply(bw_meth_list,function(x){
    makeGRangesFromDataFrame(x, seqnames.field="Chromosome", start.field="Start",  end.field="Start",
        keep.extra.columns= TRUE)
})
#add chr 
bw_meth_list <- lapply(bw_meth_list,function(x){
    seqlevelsStyle(x) <- "UCSC"
    x
})
#rename columsn
bw_meth_list <- lapply(bw_meth_list,function(x){
   colnames(mcols(x))<- c("M", "Cov", "mean_meth")
   x
})

#heaatmap
#extract matrix for plotting
ext <- 5000
mat_chip_list <- lapply(bw_chip_list, function(x){
    normalizeToMatrix(x,tss, value_column = "score" , 
    extend = ext, mean_mode = "w0", w = 10, keep = c(0, 0.99))
})
mat_meth_list <- lapply(bw_meth_list, function(x){
    normalizeToMatrix(x,tss, value_column = "mean_meth", 
    mean_mode = "absolute",
    extend = ext)
})
lapply(mat_meth_list, function(x)head(as.data.frame(x)))
#top legend
pal <- c("darkgray", "orange", "tomato3")
names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")

lgd = Legend(at = c("known",  "chimeric (novel)", "non-chimeric (novel)"), 
    title = "Transcript classification", 
    type = "lines", legend_gp = gpar(col = class_col))

#define minima and maxima
#rnaseq
common_min_rnaseq = min(unlist(mat_chip_list[grep("rnaseq", names(mat_chip_list))]))
common_max_rnaseq = max(unlist(mat_chip_list[grep("rnaseq", names(mat_chip_list))]))
common_max_ylim_rnaseq <- max(unlist(lapply(mat_chip_list[grep("rnaseq", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_rnaseq = circlize::colorRamp2(c(common_min_rnaseq, common_max_rnaseq), c("white", "forestgreen"))
#cage
common_min_cage = min(unlist(mat_chip_list[grep("cage", names(mat_chip_list))]))
common_max_cage = max(unlist(mat_chip_list[grep("cage", names(mat_chip_list))]))
common_max_ylim_cage <- max(unlist(lapply(mat_chip_list[grep("cage", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_cage = circlize::colorRamp2(c(common_min_cage, common_max_cage), c("white", "seagreen"))
#9me3
common_min_9me3 = min(unlist(mat_chip_list[grep("9me3", names(mat_chip_list))]))
common_max_9me3 = max(unlist(mat_chip_list[grep("9me3", names(mat_chip_list))]))
common_max_ylim_9me3 <- max(unlist(lapply(mat_chip_list[grep("9me3", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_9me3 = circlize::colorRamp2(c(common_min_9me3, common_max_9me3), c("white", "red2"))
#4me3
common_min_4me3 = min(unlist(mat_chip_list[grep("4me3", names(mat_chip_list))]))
common_max_4me3 = max(unlist(mat_chip_list[grep("4me3", names(mat_chip_list))]))
common_max_ylim_4me3 <- max(unlist(lapply(mat_chip_list[grep("4me3", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_4me3 = circlize::colorRamp2(c(common_min_4me3, common_max_4me3), c("white", "firebrick4"))
#H3K9ac
common_min_9ac = min(unlist(mat_chip_list[grep("9ac", names(mat_chip_list))]))
common_max_9ac = max(unlist(mat_chip_list[grep("9ac", names(mat_chip_list))]))
common_max_ylim_9ac <- max(unlist(lapply(mat_chip_list[grep("9ac", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_9ac = circlize::colorRamp2(c(common_min_9ac, common_max_9ac), c("white", "red"))
#H3K27ac
common_min_27ac = min(unlist(mat_chip_list[grep("27ac", names(mat_chip_list))]))
common_max_27ac = max(unlist(mat_chip_list[grep("27ac", names(mat_chip_list))]))
common_max_ylim_27ac <- max(unlist(lapply(mat_chip_list[grep("27ac", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_27ac = circlize::colorRamp2(c(common_min_27ac, common_max_27ac), c("white", "darkblue"))
#H3K14ac
common_min_14ac = min(unlist(mat_chip_list[grep("14ac", names(mat_chip_list))]))
common_max_14ac = max(unlist(mat_chip_list[grep("14ac", names(mat_chip_list))]))
common_max_ylim_14ac <- max(unlist(lapply(mat_chip_list[grep("14ac", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_14ac = circlize::colorRamp2(c(common_min_14ac, common_max_14ac), c("white", "cornflowerblue"))
#5ac
common_min_5ac = min(unlist(mat_chip_list[grep("5ac", names(mat_chip_list))]))
common_max_5ac = max(unlist(mat_chip_list[grep("5ac", names(mat_chip_list))]))
common_max_ylim_5ac <- max(unlist(lapply(mat_chip_list[grep("5ac", names(mat_chip_list))], function(x){
    x <- max(colMeans(x))
    x
})))
col_fun_5ac = circlize::colorRamp2(c(common_min_5ac, common_max_5ac), c("white", "#7CE3D8"))
#meth
col_fun_meth = circlize::colorRamp2(c(0, 1), c("white", "tomato3"))

#enriched heatmap creation
#top legend
ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +

    #rnaseq
    EnrichedHeatmap(mat_chip_list[["DMSO_rnaseq"]], name="DMSO_rnaseq", 
        column_title = "DMSO_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_rnaseq) +
    EnrichedHeatmap(mat_chip_list[["SB_rnaseq"]], name="SB_rnaseq", 
        column_title ="SB_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_rnaseq) +
    EnrichedHeatmap(mat_chip_list[["DAC_rnaseq"]], name="DAC_rnaseq", 
        column_title = "DAC_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_rnaseq) +
    EnrichedHeatmap(mat_chip_list[["DACSB_rnaseq"]], name="DACSB_rnaseq", 
        column_title ="DACSB_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_rnaseq) +

    #cage
    EnrichedHeatmap(mat_chip_list[["DMSO_cage"]], name="DMSO_cage", 
        column_title = "DMSO_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_cage) +
    EnrichedHeatmap(mat_chip_list[["SB_cage"]], name="SB_cage", 
        column_title ="SB_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_cage) +
    EnrichedHeatmap(mat_chip_list[["DAC_cage"]], name="DAC_cage", 
        column_title = "DAC_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_cage) +
    EnrichedHeatmap(mat_chip_list[["DACSB_cage"]], name="DACSB_cage", 
        column_title ="DACSB_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_cage) +

    #meth
    EnrichedHeatmap(mat_meth_list[["DMSO_meth"]], name="DMSO_meth", 
        column_title = "DMSO_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) +
    EnrichedHeatmap(mat_meth_list[["SB_meth"]], name="SB_meth", 
        column_title ="SB_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) +
    EnrichedHeatmap(mat_meth_list[["DAC_meth"]], name="DAC_meth", 
        column_title = "DAC_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) +
    EnrichedHeatmap(mat_meth_list[["DACSB_meth"]], name="DACSB_meth", 
        column_title ="DACSB_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) +

    #9me3
    EnrichedHeatmap(mat_chip_list[["DMSO_9me3"]], name="DMSO_9me3", 
        column_title = "DMSO_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9me3) +
    EnrichedHeatmap(mat_chip_list[["SB_9me3"]], name="SB_9me3", 
        column_title ="SB_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9me3) +
    EnrichedHeatmap(mat_chip_list[["DAC_9me3"]], name="DAC_9me3", 
        column_title = "DAC_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9me3) +
    EnrichedHeatmap(mat_chip_list[["DACSB_9me3"]], name="DACSB_9me3", 
        column_title ="DACSB_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9me3) +

    #4me3
    EnrichedHeatmap(mat_chip_list[["DMSO_4me3"]], name="DMSO_4me3", 
        column_title = "DMSO_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_4me3) +
    EnrichedHeatmap(mat_chip_list[["SB_4me3"]], name="SB_4me3", 
        column_title ="SB_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_4me3) +
    EnrichedHeatmap(mat_chip_list[["DAC_4me3"]], name="DAC_4me3", 
        column_title = "DAC_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_4me3) +
    EnrichedHeatmap(mat_chip_list[["DACSB_4me3"]], name="DACSB_4me3", 
        column_title ="DACSB_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_4me3) +

    #9ac
    EnrichedHeatmap(mat_chip_list[["DMSO_9ac"]], name="DMSO_9ac", 
        column_title = "DMSO_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9ac) +
    EnrichedHeatmap(mat_chip_list[["SB_9ac"]], name="SB_9ac", 
        column_title ="SB_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_9ac"]], name="DAC_9ac", 
        column_title = "DAC_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_9ac"]], name="DACSB_9ac", 
        column_title ="DACSB_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9ac) +

    #27ac
    EnrichedHeatmap(mat_chip_list[["DMSO_27ac"]], name="DMSO_27ac", 
        column_title = "DMSO_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_27ac) +
    EnrichedHeatmap(mat_chip_list[["SB_27ac"]], name="SB_27ac", 
        column_title ="SB_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_27ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_27ac"]], name="DAC_27ac", 
        column_title = "DAC_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_27ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_27ac"]], name="DACSB_27ac", 
        column_title ="DACSB_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_27ac) +
    
    #14ac
    EnrichedHeatmap(mat_chip_list[["DMSO_14ac"]], name="DMSO_14ac", 
        column_title = "DMSO_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_14ac) +
    EnrichedHeatmap(mat_chip_list[["SB_14ac"]], name="SB_14ac", 
        column_title ="SB_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_14ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_14ac"]], name="DAC_14ac", 
        column_title = "DAC_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_14ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_14ac"]], name="DACSB_14ac", 
        column_title ="DACSB_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_14ac) +

    #5ac
    EnrichedHeatmap(mat_chip_list[["DMSO_5ac"]], name="DMSO_5ac", 
        column_title = "DMSO_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) +
    EnrichedHeatmap(mat_chip_list[["SB_5ac"]], name="SB_5ac", 
        column_title ="SB_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_5ac"]], name="DAC_5ac", 
        column_title = "DAC_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_5ac"]], name="DACSB_5ac", 
        column_title ="DACSB_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) 


#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_all.pdf"), width=38, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),9)), "mm"))
dev.off() 





#only rnaseq
ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +

    #rnaseq
    EnrichedHeatmap(mat_chip_list[["DMSO_rnaseq"]], name="DMSO_rnaseq", 
        column_title = "DMSO_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_rnaseq) +
    EnrichedHeatmap(mat_chip_list[["SB_rnaseq"]], name="SB_rnaseq", 
        column_title ="SB_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_rnaseq) +
    EnrichedHeatmap(mat_chip_list[["DAC_rnaseq"]], name="DAC_rnaseq", 
        column_title = "DAC_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_rnaseq) +
    EnrichedHeatmap(mat_chip_list[["DACSB_rnaseq"]], name="DACSB_rnaseq", 
        column_title ="DACSB_rnaseq",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_rnaseq, accuracy=ifelse(common_max_ylim_rnaseq<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_rnaseq) 

#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_rnaseq.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 



#only cage
ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #cage
    EnrichedHeatmap(mat_chip_list[["DMSO_cage"]], name="DMSO_cage", 
        column_title = "DMSO_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_cage) +
    EnrichedHeatmap(mat_chip_list[["SB_cage"]], name="SB_cage", 
        column_title ="SB_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_cage) +
    EnrichedHeatmap(mat_chip_list[["DAC_cage"]], name="DAC_cage", 
        column_title = "DAC_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_cage) +
    EnrichedHeatmap(mat_chip_list[["DACSB_cage"]], name="DACSB_cage", 
        column_title ="DACSB_cage",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_cage, accuracy=ifelse(common_max_ylim_cage<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_cage) 
#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_cage.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 




#only meth
ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #meth
    EnrichedHeatmap(mat_meth_list[["DMSO_meth"]], name="DMSO_meth", 
        column_title = "DMSO_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) +
    EnrichedHeatmap(mat_meth_list[["SB_meth"]], name="SB_meth", 
        column_title ="SB_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) +
    EnrichedHeatmap(mat_meth_list[["DAC_meth"]], name="DAC_meth", 
        column_title = "DAC_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) +
    EnrichedHeatmap(mat_meth_list[["DACSB_meth"]], name="DACSB_meth", 
        column_title ="DACSB_meth",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,1))), col = col_fun_meth) 
#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_meth.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 


#only 9me3

ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #9me3
    EnrichedHeatmap(mat_chip_list[["DMSO_9me3"]], name="DMSO_9me3", 
        column_title = "DMSO_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9me3) +
    EnrichedHeatmap(mat_chip_list[["SB_9me3"]], name="SB_9me3", 
        column_title ="SB_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9me3) +
    EnrichedHeatmap(mat_chip_list[["DAC_9me3"]], name="DAC_9me3", 
        column_title = "DAC_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9me3) +
    EnrichedHeatmap(mat_chip_list[["DACSB_9me3"]], name="DACSB_9me3", 
        column_title ="DACSB_9me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9me3, accuracy=ifelse(common_max_ylim_9me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9me3) 
#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_9me3.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 


#only 4me3
ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #4me3
    EnrichedHeatmap(mat_chip_list[["DMSO_4me3"]], name="DMSO_4me3", 
        column_title = "DMSO_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_4me3) +
    EnrichedHeatmap(mat_chip_list[["SB_4me3"]], name="SB_4me3", 
        column_title ="SB_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_4me3) +
    EnrichedHeatmap(mat_chip_list[["DAC_4me3"]], name="DAC_4me3", 
        column_title = "DAC_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_4me3) +
    EnrichedHeatmap(mat_chip_list[["DACSB_4me3"]], name="DACSB_4me3", 
        column_title ="DACSB_4me3",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_4me3, accuracy=ifelse(common_max_ylim_4me3<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_4me3) 
#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_4me3.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 


#9ac

ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #9ac
    EnrichedHeatmap(mat_chip_list[["DMSO_9ac"]], name="DMSO_9ac", 
        column_title = "DMSO_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9ac) +
    EnrichedHeatmap(mat_chip_list[["SB_9ac"]], name="SB_9ac", 
        column_title ="SB_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_9ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_9ac"]], name="DAC_9ac", 
        column_title = "DAC_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_9ac"]], name="DACSB_9ac", 
        column_title ="DACSB_9ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_9ac, accuracy=ifelse(common_max_ylim_9ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_9ac) 
#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_9ac.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 


#only 27 ac
ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #27ac
    EnrichedHeatmap(mat_chip_list[["DMSO_27ac"]], name="DMSO_27ac", 
        column_title = "DMSO_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_27ac) +
    EnrichedHeatmap(mat_chip_list[["SB_27ac"]], name="SB_27ac", 
        column_title ="SB_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_27ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_27ac"]], name="DAC_27ac", 
        column_title = "DAC_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_27ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_27ac"]], name="DACSB_27ac", 
        column_title ="DACSB_27ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_27ac, accuracy=ifelse(common_max_ylim_27ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_27ac) 

#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_27ac.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 


#only 10ac
ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #14ac
    EnrichedHeatmap(mat_chip_list[["DMSO_14ac"]], name="DMSO_14ac", 
        column_title = "DMSO_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_14ac) +
    EnrichedHeatmap(mat_chip_list[["SB_14ac"]], name="SB_14ac", 
        column_title ="SB_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_14ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_14ac"]], name="DAC_14ac", 
        column_title = "DAC_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_14ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_14ac"]], name="DACSB_14ac", 
        column_title ="DACSB_14ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_14ac, accuracy=ifelse(common_max_ylim_14ac<10, 1,10)*1.25,
        f = ceiling)))), 
        col =  col_fun_14ac) 
#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_14ac.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 


#only 5ac

ht_list <-  Heatmap( as.vector(tss$class_code_simple), col = class_col, name = "Transcript classification",
    show_row_names = FALSE, width = unit(3, "mm")) +
    #5ac
    EnrichedHeatmap(mat_chip_list[["DMSO_5ac"]], name="DMSO_5ac", 
        column_title = "DMSO_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) +
    EnrichedHeatmap(mat_chip_list[["SB_5ac"]], name="SB_5ac", 
        column_title ="SB_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) +
    EnrichedHeatmap(mat_chip_list[["DAC_5ac"]], name="DAC_5ac", 
        column_title = "DAC_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) +
    EnrichedHeatmap(mat_chip_list[["DACSB_5ac"]], name="DACSB_5ac", 
        column_title ="DACSB_5ac",
        clustering_distance_rows = dist_by_closeness,
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = class_col), 
        ylim=c(0,plyr::round_any(common_max_ylim_5ac, accuracy=ifelse(common_max_ylim_5ac<10, 1,10)*1.25, 
        f = ceiling)))), 
        col =  col_fun_5ac) 


#plot heatmap
#test
dir.create(file.path(PostDE.dir, "enrichedHeat"),recursive=TRUE)
pdf(file.path(PostDE.dir, "enrichedHeat", "enrichedHeat_5ac.pdf"), width=10, height=7)
draw(ht_list, split = factor( as.vector(tss$class_code_simple),c( "chimeric (novel)", "non-chimeric (novel)","known")),cluster_rows = F,
    ht_gap = unit(c(2, rep(c(6,6,6,6),1)), "mm"))
dev.off() 
