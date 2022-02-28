#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)qq
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)


#load H1299 data
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
tpm_mean_h1299 <- readRDS(file.path(results.dir, "tpm_meanTreatment_from_counts.rds"))

#Set specifications
alpha <- 0.01 #set FDR cutoff
cutoff <- alpha
lfc <- 2##set logfold2 cutoff
l2fc <- lfc

#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id 

#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]

#merge info
tpm_mean_h1299 <- as.data.frame(tpm_mean_h1299)
tpm_mean_h1299$transcript_id <- rownames(tpm_mean_h1299)
tpm_mean_h1299_merged <- dplyr::left_join(tpm_mean_h1299, anno_classi[, c("transcript_id", "class_code_simple")])
tpm_mean_h1299_merged$transcript_id <- NULL
tpm_mean_h1299_merged_melt <- reshape2::melt(tpm_mean_h1299_merged)
tpm_mean_h1299_merged_melt$class_code_simple<- factor(tpm_mean_h1299_merged_melt$class_code_simple, 
    levels=rev(c("non-chimeric (novel)", "chimeric (novel)" ,"known")))
pal <- c("darkgray", "orange", "tomato3")
font_size <- c(16,"plain", "black")
names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
tpm_allTranscripts <- ggboxplot(tpm_mean_h1299_merged_melt, x="variable", y="value",
    fill="class_code_simple",  palette=class_col,lab.pos="out",#color="Var1",
    ylab="TPM",  order=rev(c("DMSO", "SB939", "DAC", "DACandSB939")),outlier.shape=NA,
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    ) +
    rremove("legend.title") +
    rremove("ylab") +
    coord_flip() +
    ylim(c(0,10))
pdf(file.path(PreDE.dir, "transcript_classification_TPM.pdf"), height=5, width=5 )
print(tpm_allTranscripts)
dev.off()