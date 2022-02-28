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
repeats <- frea d("/omics/groups/OE0219/internal/repguide/data/repeats/hg19_repeats.txt")

#rnaseq
bw_rna <- list("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DMSO.normalized.bigWig","/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/SB939.normalized.bigWig",
    "/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DAC.normalized.bigWig","/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DACSB939.normalized.bigWig")

#cage
bw_cage <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DMSO.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/SB939.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DAC.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DACSB.bigWig")

#H3K9me3
bw_9me3 <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150368_DMSO_H3K9me3_ATCACG_L002_R1_complete_filtered.fastq.gz.normalized.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150384_SB939_H3K9me3_TTAGGC_L002_R1_complete_filtered.fastq.gz.normalized.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150352_DAC_H3K9me3_CGATGT_L002_R1_complete_filtered.fastq.gz.normalized.bigWig",  "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150360_DACSB_H3K9me3_ACTTGA_L002_R1_complete_filtered.fastq.gz.normalized.bigWig")

#H3K4me3
bw_4me3 <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150366_DMSO_H3K4me3_ATCACG_L001_R1_complete_filtered.fastq.gz.normalized.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150382_SB939_H3K4me3_ACTTGA_L001_R1_complete_filtered.fastq.gz.normalized.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150350_DAC_H3K4me3_GCCAAT_L001_R1_complete_filtered.fastq.gz.normalized.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150358_DACSB_H3K4me3_CTTGTA_L001_R1_complete_filtered.fastq.gz.normalized.bigWig")

#H3K9ac
bw_9ac <- list( "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150367_DMSO_H3K9ac_CGATGT_L005_R1_complete_filtered.fastq.gz.normalized.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150383_SB939_H3K9ac_TGACCA_L005_R1_complete_filtered.fastq.gz.normalized.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150351_DAC_H3K9ac_TTAGGC_L005_R1_complete_filtered.fastq.gz.normalized.bigWig",  "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150359_DACSB_H3K9ac_GATCAG_L005_R1_complete_filtered.fastq.gz.normalized.bigWig")
#H3K27ac
bw_27ac <- list( "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150362_DMSO_H3K27ac_ATCACG_L008_R1_complete_filtered.fastq.gz.normalized.bigWig",  "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150378_SB939_H3K27ac_ACTTGA_L008_R1_complete_filtered.fastq.gz.normalized.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150346_DAC_H3K27ac_GCCAAT_L008_R1_complete_filtered.fastq.gz.normalized.bigWig",  "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150354_DACSB_H3K27ac_CTTGTA_L008_R1_complete_filtered.fastq.gz.normalized.bigWig")

#load chip data part 1
input.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/set_1"
sample_anno_chip <- fread(file.path(input.dir, "samples.csv"))
sample_anno_chip <- sample_anno_chip[sample_anno_chip$Modification != "Input",]
sample_anno_chip$normalizedBigWig <- file.path(input.dir, paste0(sample_anno_chip$"Sample ID",".normalized.bigWig"))
sample_anno_chip_split <- split(sample_anno_chip,sample_anno_chip$Modification)
sample_anno_chip_split <- lapply(sample_anno_chip_split, function(x){
    y <- x$normalizedBigWig
    names(y) <-x$Treatment
    y <- y[c("DMSO", "SB", "DAC", "DAC+SB")]
}
)

#load chip data part 2
#input.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/set_2"
#files <- list.files(input.dir, pattern="normalized.bigWig", full.names=TRUE)
# sample_anno_chip_split_2 <- list(H3K27ac=files[grep("H3K27ac", files)][c(3,1,5,2)],
#     H3K4me1=files[grep("H3K4me1", files)][c(3,1,5,2)],  H3K4me3=files[grep("H3K4me3", files)][c(3,1,5,2)],  
#     H3K9ac=files[grep("H3K9ac", files)][c(3,1,5,2)])
#sample_anno_chip_split_2 <- list(H3K4me1=files[grep("H3K4me1", files)][c(3,1,5,2)])

# #load chip data part 3
# input.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/set_3"
# files <- list.files(input.dir, pattern="bigWig", full.names=TRUE)
# sample_anno_chip_split_2 <- list(H3K27ac=files[grep("H3K27ac", files)][c(3,1,5,2)],
#     H3K4me1=files[grep("H3K4me1", files)][c(3,1,5,2)],  H3K4me3=files[grep("H3K4me3", files)][c(3,1,5,2)],  
#     H3K9ac=files[grep("H3K9ac", files)][c(3,1,5,2)])
    
#combine lists
#bw_list <- c(list(bw_rna, bw_cage), sample_anno_chip_split, sample_anno_chip_split_2)
#names(bw_list) <- c("RNAseq","CAGE", names(sample_anno_chip_split), names(sample_anno_chip_split_2))
bw_list <- c(list(bw_rna, bw_cage, bw_9me3,bw_4me3, bw_9ac , bw_27ac), sample_anno_chip_split)
names(bw_list) <- c("RNAseq","CAGE","H3K9me3","H3K4me3", "H3K9ac","H3K27ac",   names(sample_anno_chip_split))

#prepare lists and sample sheets
sample_anno <- data.frame(sample_name = factor(c("DMSO", "SB", "DAC", "DACSB"),levels=c("DMSO", "SB", "DAC", "DACSB")), treatment=factor(c("DMSO", "SB", "DAC", "DACSB"),levels=c("DMSO", "SB", "DAC", "DACSB")), bigwig=NA)

#prepare lists and sample sheets 
col <- c("#FC4E07","gray","#00AFBB", "#E7B800")
names(col) <-c("DMSO", "SB","DAC", "DACSB")
#get regions of interest 
anno_classi$transcript_id <- anno_classi$qry_id
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
DEG_results_list_anno_sub <- lapply(DEG_results_list, function(x){
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    x <- split(x, x$class_code_simple)
    x
})
DEG_results_list_anno_sub_sub <- DEG_results_list_anno_sub$DACandSB939_vs_DMSO

#loop over different sets of interest and plot them
dir.create(file.path(PostDE.dir, "profilePlots", "data"), recursive=TRUE)
dir.create(file.path(PostDE.dir, "profilePlots", "normalized"))

anno_transcript <- anno[anno$type=="transcript",]

for(ext in c( 5000, 500, 1500)){
    for(i in names(DEG_results_list_anno_sub_sub)){
        trans <- DEG_results_list_anno_sub_sub[[i]]$transcript_id
        tss <- resize(anno_transcript[anno_transcript$transcript_id %in% trans, ],1)
        mcols(tss) <- NULL
        export.bed(tss, file.path(PostDE.dir, "profilePlots", "data", "tss_oi.bed"))

        for(j in names(bw_list)[3]){
            sample_anno$bigwig <- unlist(bw_list[[j]])
            #create column data
            col_data <- read_coldata(coldata=sample_anno, files_idx=3, sample_idx=1)
            #create matrix list
            regions_matrix <- extract_matrix(coldata=col_data, 
                bed=  file.path(PostDE.dir, "profilePlots", "data", "tss_oi.bed"), 
                genome=hg19, up=ext, down=ext, binSize=10, 
                nthreads=4)
        #draw profile plot 
            pdf(file.path(PostDE.dir, "profilePlots","normalized", paste(i, "_", j,"_profilePlot_ext", ext,".pdf")), height=3,5, width=3.5)
                print(plot_profile(mat_list=regions_matrix, condition="treatment", 
                collapse_reps=TRUE,summarizeBy="mean", ci=TRUE, color=col))
            dev.off()
            print(paste0(j))
        }
        print(paste0(i))
    }
        print(paste0(ext))
}














# #heaatmap
# #extract matrix for plotting
# ext <- 500
# mat_list <- lapply(bw_list, function(x){
#     normalizeToMatrix(x,repeats_LTR5all, value_column = "score" , 
#     extend = ext, mean_mode = "w0", w = 10, keep = c(0, 0.99))
# })
# #common color annotation
# common_min = min(unlist(mat_list))
# common_max = max(unlist(mat_list))
# col_fun = circlize::colorRamp2(c(common_min, common_max), c("white", "red"))
# #top legend
# pal <- c("darkgray", "orange", "tomato3")
# names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")

# lgd = Legend(at = c("LTR5_Hs", "LTR5B", "LTR5A"), title = "Repeat names", 
#     type = "lines", legend_gp = gpar(col = 2:4))


# #enriched heatmap creation
# ht_list <-  Heatmap( as.vector(repeats_LTR5all$repName), col = structure(2:4, names =c("LTR5_Hs", "LTR5B", "LTR5A")), name = "Repeat names",
#     show_row_names = FALSE, width = unit(3, "mm")) +
    
#     EnrichedHeatmap(mat_list[[1]], name=names(mat_list)[[1]], column_title = "Fuentes et al. Sp gRNAs\noriginal data",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +
    
#     EnrichedHeatmap(mat_list[[2]], name=names(mat_list)[[2]], column_title = "Non-targeting Sp gRNAs\noriginal data",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[3]], name=names(mat_list)[[3]], column_title = "Fuentes et al. Sa gRNAs\noriginal data",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[4]], name=names(mat_list)[[4]], column_title = "Fuentes et al. gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[5]], name=names(mat_list)[[5]], column_title = "RepGuide loose gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means), 
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[6]], name=names(mat_list)[[6]], column_title = "RepGuide strict gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))),
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[7]], name=names(mat_list)[[7]], column_title = "Non-targeting gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))),  
#     col =  col_fun) 

# #plot heatmap
# dir.create(file.path(analysis.dir, "bw_heatmaps", "merged_bw"),recursive=TRUE)
# pdf(file.path(analysis.dir, "bw_heatmaps","merged_bw",paste0("enrichedHeat_bwMergedLTR5_repNameStrat_noCluster_bin10_ext", ext,"_adjFuentes.pdf")), width=20, height=7)
# draw(ht_list, split =  as.vector(repeats_LTR5all$repName),cluster_rows = F,
#     ht_gap = unit(c(2, 6, 6,6, 6, 6,6), "mm"))
# dev.off() 