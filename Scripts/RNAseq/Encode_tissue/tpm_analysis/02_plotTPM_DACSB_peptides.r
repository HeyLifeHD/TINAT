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
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"

#Read in Data
vst_h1299 <- readRDS(file.path(results_h1299.dir, "vst.rds"))
dds_h1299 <-readRDS(file.path(results_h1299.dir, "dds.rds"))
DEG_results_list_h1299<- readRDS(file.path(PostDE_h1299.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
tpm_mean_h1299 <- readRDS(file.path(results_h1299.dir, "tpm_meanTreatment_from_counts.rds"))
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))

#Set specifications
alpha <- 0.01 #set FDR cutoff
cutoff <- alpha
lfc <- 2##set logfold2 cutoff
l2fc <- lfc

#get candidate transcripts
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
peptides_new_sub <- peptides_new_ORF[peptides_new_ORF$DAC_SB==100 & peptides_new_ORF$DMSO ==0,]

tpm_mean_combined <- readRDS(file.path(results.dir, "tpm_mean_combined_from_counts.rds"))
#combine annotation
col <- randomcoloR::distinctColorPalette(length(unique(colData(dds)$tissue)))
col <- c(col, c("#00AFBB", "#E7B800", "#FC4E07","gray"))
names(col) <- c(levels(colData(dds)$tissue), levels(colData(dds_h1299)$treatment))
anno_colors <- list(tissue=col)
annovst_tissue <- data.frame(row.names= c(levels(anno$tissue),levels(colData(dds_h1299)$treatment)), 
tissue= c(levels(anno$tissue),levels(colData(dds_h1299)$treatment)))


#plot them
matvst <- log10(tpm_mean_combined[peptides_new_sub$transcript_id    ,]+1)
pdf(file.path(PreDE.dir,"TPM_log10_meanCombined_Heatmap_candidateTranscripts_absolute_noScale.pdf"),height= 7)
pheatmap(matvst, scale="none", show_colnames=T,color= c(rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()
matvst <- tpm_mean_combined[peptides_new_sub$transcript_id    ,]
pdf(file.path(PreDE.dir,"TPM_meanCombined_Heatmap_candidateTranscripts_absolute_rowScale.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=T,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst_tissue), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()