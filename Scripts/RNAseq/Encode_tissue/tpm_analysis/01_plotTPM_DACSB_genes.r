#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)qq
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
tpm <- readRDS(file.path(results.dir, "tpm_from_counts.rds"))
tpm_mean <- readRDS(file.path(results.dir, "tpm_meanTissue_from_counts.rds"))


#load H1299 data
#Directories
base_h1299.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results_h1299.dir <- file.path(base_h1299.dir, "results")
results_h1299.dir<- file.path(base_results_h1299.dir , "tables")
PreDE_h1299.dir <- file.path(base_results_h1299.dir,"PreDE")
PostDE_h1299.dir <- file.path(base_results_h1299.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

#Read in Data
vst_h1299 <- readRDS(file.path(results_h1299.dir, "vst.rds"))
dds_h1299 <-readRDS(file.path(results_h1299.dir, "dds.rds"))
DEG_results_list_h1299<- readRDS(file.path(PostDE_h1299.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/cinternal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
tpm_mean_h1299 <- readRDS(file.path(results_h1299.dir, "tpm_meanTreatment_from_counts.rds"))

#Set specifications
alpha <- 0.01 #set FDR cutoff
cutoff <- alpha
lfc <- 2##set logfold2 cutoff
l2fc <- lfc

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

#plotting parameters
#sample level
col <- randomcoloR::distinctColorPalette(length(colData(dds)$tissue))
names(col) <- levels(colData(dds)$tissue)
anno_colors <- list(tissue=col)
anno <- colData(vst)
annovst <- anno[,c("tissue" ), drop=FALSE]
#tissue level

annovst_tissue <- data.frame(row.names= unique(anno$tissue), tissue= unique(anno$tissue))
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
matvst <- tpm[genes2plot$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_Heatmap_all_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()
matvst <- tpm_mean[genes2plot$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_meanTissue_Heatmap_all_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

##get upregulated -novel DEGs
genes2plot_novel <- lapply(DEG_results_list_sub_up, function(x){
  x <- x[x$class_code_simple != "known",]
  x <- x$transcript_id
  x
})
length(genes2plot_novel$DACandSB939_vs_DMSO)
#plot them
matvst <- tpm[genes2plot_novel$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_Heatmap_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()
matvst <- tpm_mean[genes2plot_novel$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_meanTissue_Heatmap_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
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
matvst <- tpm[genes2plot_novel_LTR12$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_Heatmap_LTR12_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()
matvst <- tpm_mean[genes2plot_novel_LTR12$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_meanTissue_Heatmap_LTR12_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()



#combined visualization
#combine tpm
tpm_mean <- as.data.frame(tpm_mean)
tpm_mean$transcript_id <- rownames(tpm_mean)
tpm_mean_h1299 <- as.data.frame(tpm_mean_h1299)
tpm_mean_h1299$transcript_id <- rownames(tpm_mean_h1299)
tpm_mean_combined <- dplyr::inner_join(tpm_mean, tpm_mean_h1299, by="transcript_id")
rownames(tpm_mean_combined)<- tpm_mean_combined$transcript_id
tpm_mean_combined$transcript_id <- NULL
saveRDS(tpm_mean_combined, file.path(results.dir, "tpm_mean_combined_from_counts.rds"))
data.table::fwrite(tpm_mean_combined, file.path(results.dir, "tpm_mean_combined_from_counts.txt"), row.names=TRUE)
#combine annotation
col <- randomcoloR::distinctColorPalette(length(unique(colData(dds)$tissue)))
col <- c(col, treat_col)
names(col) <- c(levels(colData(dds)$tissue), names(treat_col))
anno_colors <- list(tissue=col)
annovst_tissue <- data.frame(row.names= c(levels(anno$tissue),levels(colData(dds_h1299)$treatment)), 
tissue= c(levels(anno$tissue),levels(colData(dds_h1299)$treatment)))

#get upregulated DEGs
genes2plot <- lapply(DEG_results_list_sub_up, function(x){
  x <- x$transcript_id
  x
})
#plot them
matvst <- log10(tpm_mean_combined[genes2plot$DACandSB939_vs_DMSO,]+1)
pdf(file.path(PreDE.dir,"TPM_log10_meanCombined_Heatmap_all_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()
matvst <- tpm_mean_combined[genes2plot$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_meanCombined_Heatmap_all_up_DACSB_absolute_rowScale.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

##get upregulated -novel DEGs
genes2plot_novel <- lapply(DEG_results_list_sub_up, function(x){
  x <- x[x$class_code_simple != "known",]
  x <- x$transcript_id
  x
})
length(genes2plot_novel$DACandSB939_vs_DMSO)
#plot them
matvst <-  log10(tpm_mean_combined[genes2plot_novel$DACandSB939_vs_DMSO,]+1)
pdf(file.path(PreDE.dir,"TPM_log10_meanCombined_Heatmap_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()
matvst <- tpm_mean_combined[genes2plot_novel$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_meanCombined_Heatmap_novel_up_DACSB_absolute_rowScale.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()


#as boxplot
matvst <- tpm_mean_combined[genes2plot_novel$DACandSB939_vs_DMSO,]
mat_plot <- matvst
mat_plot <- reshape2::melt(mat_plot)
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_novel_up_DACSB_absolute_withoutOutliers.pdf"),height= 10)
ggboxplot(mat_plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, differentially induced transcripts", fill="variable", palette=anno_colors$tissue, outlier.shape=NA)+
    rremove("ylab")+ylim(c(0,20))+rremove("legend")
dev.off()
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_novel_up_DACSB_absolute.pdf"),height= 10)
ggboxplot(mat_plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel,differentially induced transcripts", fill="variable", palette=anno_colors$tissue)+
    rremove("ylab")+rremove("legend")
dev.off()

#combine with housekeeping list
hkt <- data.table::fread(file.path(results.dir,"MostStable.csv" ))
#get our transcript anno
anno_classi_known <- anno_classi[anno_classi$class_code =="=", ]
hkt_own <- anno_classi_known$transcript_id[sapply(strsplit(as.character(anno_classi_known$ref_id), ".", fixed=TRUE), "[",1) %in% hkt$"Ensembl ID" ]
length(hkt_own)
house_keeper <- tpm_mean_combined[ hkt_own,]
house_keeper <- reshape2::melt(house_keeper)
mat_plot$class <- "novel, upregulated, and LTR12-derived DEGs"
house_keeper$class <- "housekeping transcripts"
 plot <- rbind(mat_plot, house_keeper)

pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_novel_up_DACSB_absolute_withoutOutliers_withHKT.pdf"),height= 10)
ggboxplot(plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, upregulated, and LTR12-derived DEGs", fill="class", outlier.shape=NA)+
    rremove("ylab")+ylim(c(0,20))+rremove("legend")
dev.off()
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_novel_up_DACSB_absolute_withHKT.pdf"),height= 10)
ggboxplot(plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, differentially induced", fill="class")+
    rremove("ylab")+rremove("legend")
dev.off()

#only house keeping
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot__novel_up_DACSB_absolute_withoutOutliers_onlyHKT.pdf"),height= 10)
ggboxplot(house_keeper, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of housekeeper transcripts", fill="variable", outlier.shape=NA, palette=anno_colors$tissue)+
    rremove("ylab")+ylim(c(0,20))+rremove("legend")
dev.off()
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_novel_up_DACSB_absolute_onlyHKT.pdf"),height= 10)
ggboxplot(house_keeper, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of housekeeper transcripts", fill="variable", palette=anno_colors$tissue)+
    rremove("ylab")+rremove("legend")
dev.off()




#genes upregulated, novel and ltr12 derived
genes2plot_novel_LTR12 <- lapply(DEG_results_list_sub_up, function(x){
  x <- x[x$class_code_simple != "known",]
  x <- x[which(x$dist_nearest_LTR12repeat ==0),]
  x <- x$transcript_id
  x
})
length(genes2plot_novel_LTR12$DACandSB939_vs_DMSO)
#plot them
matvst <-  log10(tpm_mean_combined[genes2plot_novel_LTR12$DACandSB939_vs_DMSO,]+1)
pdf(file.path(PreDE.dir,"TPM_log10_meanCombined_Heatmap_LTR12_novel_up_DACSB_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=F,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=F,annotation_colors=anno_colors, fontsize_row=5)
dev.off()
matvst <- tpm_mean_combined[genes2plot_novel_LTR12$DACandSB939_vs_DMSO,]
pdf(file.path(PreDE.dir,"TPM_meanCombined_Heatmap_LTR12_novel_up_DACSB_absolute_rowScale.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=T,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=F,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

#as boxplot
matvst <- tpm_mean_combined[genes2plot_novel_LTR12$DACandSB939_vs_DMSO,]
mat_plot <- matvst
mat_plot <- reshape2::melt(mat_plot)
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_LTR12_novel_up_DACSB_absolute_withoutOutliers.pdf"),height= 10)
ggboxplot(mat_plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, upregulated, and LTR12-derived DEGs", fill="variable", palette=anno_colors$tissue, outlier.shape=NA)+
    rremove("ylab")+ylim(c(0,20))+rremove("legend")
dev.off()
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_LTR12_novel_up_DACSB_absolute.pdf"),height= 10)
ggboxplot(mat_plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, upregulated, and LTR12-derived DEGs", fill="variable", palette=anno_colors$tissue)+
    rremove("ylab")+rremove("legend")
dev.off()

#combine with housekeeping list
hkt <- data.table::fread(file.path(results.dir,"MostStable.csv" ))
#get our transcript anno
anno_classi_known <- anno_classi[anno_classi$class_code =="=", ]
hkt_own <- anno_classi_known$transcript_id[sapply(strsplit(as.character(anno_classi_known$ref_id), ".", fixed=TRUE), "[",1) %in% hkt$"Ensembl ID" ]
length(hkt_own)
house_keeper <- tpm_mean_combined[ hkt_own,]
house_keeper <- reshape2::melt(house_keeper)
mat_plot$class <- "novel, upregulated, and LTR12-derived DEGs"
house_keeper$class <- "housekeping transcripts"
 plot <- rbind(mat_plot, house_keeper)

pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_LTR12_novel_up_DACSB_absolute_withoutOutliers_withHKT.pdf"),height= 10)
ggboxplot(plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, upregulated, and LTR12-derived DEGs", fill="class", outlier.shape=NA)+
    rremove("ylab")+ylim(c(0,20))+rremove("legend")
dev.off()
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_LTR12_novel_up_DACSB_absolute_withHKT.pdf"),height= 10)
ggboxplot(plot, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, upregulated, and LTR12-derived DEGs", fill="class")+
    rremove("ylab")+rremove("legend")
dev.off()

#only house keeping
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_LTR12_novel_up_DACSB_absolute_withoutOutliers_onlyHKT.pdf"),height= 10)
ggboxplot(house_keeper, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, upregulated, and LTR12-derived DEGs", fill="variable", outlier.shape=NA, palette=anno_colors$tissue)+
    rremove("ylab")+ylim(c(0,20))+rremove("legend")
dev.off()
pdf(file.path(PreDE.dir,"TPM_meanCombined_boxplot_LTR12_novel_up_DACSB_absolute_onlyHKT.pdf"),height= 10)
ggboxplot(house_keeper, x="variable", y="value",  order=rev(colnames(matvst[,order(colMedians(as.matrix(matvst)), decreasing=TRUE)])), 
    rotate=TRUE, ylab="TPM of novel, upregulated, and LTR12-derived DEGs", fill="variable", palette=anno_colors$tissue)+
    rremove("ylab")+rremove("legend")
dev.off()