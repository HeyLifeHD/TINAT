#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)

#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/Encode/220128_deNovo_quantification_H1299ref_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#DESeq2 Analysis
dds <-readRDS(file.path(results.dir, "dds.rds"))
vst <- readRDS(file.path(results.dir,"vst.rds"))

#load H1299 data
#Directories
base_h1299.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results_h1299.dir <- file.path(base_h1299.dir, "results")
results_h1299.dir<- file.path(base_results.dir , "tables")
PreDE_h1299.dir <- file.path(base_results_h1299.dir,"PreDE")
PostDE_h1299.dir <- file.path(base_results_h1299.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

#Read in Data
vst_h1299 <- readRDS(file.path(results_h1299.dir, "vst.rds"))
dds_h1299 <-readRDS(file.path(results_h1299.dir, "dds.rds"))
DEG_results_list_h1299<- readRDS(file.path(PostDE_h1299.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))

#Set specifications
alpha <- 0.01 #set FDR cutoff
cutoff <- alpha
lfc <- 2##set logfold2 cutoff
l2fc <- lfc

#plotting parameters
col <- randomcoloR::distinctColorPalette(length(colData(dds)$tissue))
names(col) <- levels(colData(dds)$tissue)
anno_colors <- list(tissue=col)
anno <- colData(vst)
annovst <- anno[,c("tissue" ), drop=FALSE]

#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id 

#get transcript and repeat annotation
anno_classi$transcript_id <- anno_classi$qry_id 
DEG_results_list_sub <- DEG_results_list_h1299[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
DEG_results_list_sub_up <- lapply(DEG_results_list_sub, function(x){
    x<- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x$ERE_TSS_anno <- NA
    x$ERE_TSS_anno <- ifelse(x$dist_nearest_repeat == 0, x$nearest_repeat_repClass, "no ERE")
    x$ERE_TSS_anno <- ifelse(x$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), x$ERE_TSS_anno , "other")
    x
})

#get upregulated DEGs
genes2plot <- lapply(DEG_results_list_sub_up, function(x){
  x <- x$transcript_id
  x
})
#plot them
matvst <- assay(vst)[genes2plot$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"Heatmap_all_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

##get upregulated -novel DEGs
genes2plot_novel <- lapply(DEG_results_list_sub_up, function(x){
  x <- x[x$class_code_simple != "known",]
  x <- x$transcript_id
  x
})
length(genes2plot_novel$DACandSB939_vs_DMSO)
#plot them
matvst <- assay(vst)[genes2plot_novel$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"Heatmap_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

#genes upregulated, novel and ltr12 derived
##get upregulated -novel DEGs
genes2plot_novel_LTR12 <- lapply(DEG_results_list_sub_up, function(x){
  x <- x[x$class_code_simple != "known",]
  x <- x[which(x$dist_nearest_LTR12repeat ==0),]
  x <- x$transcript_id
  x
})
length(genes2plot_novel_LTR12$DACandSB939_vs_DMSO)
#plot them
matvst <- assay(vst)[genes2plot_novel_LTR12$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"Heatmap_LTR12_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()