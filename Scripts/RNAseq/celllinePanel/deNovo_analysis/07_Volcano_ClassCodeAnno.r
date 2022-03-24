#Investigate overlap of DEGs with transcript anno
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(knitr)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)
library(clusterProfiler)
library(readxl)
library(rtracklayer)
library(LOLA)
library(ChIPpeakAnno)
library(rtracklayer)
library(ggpubr)
library(LOLA)
library(GenomicFeatures)
library(UpSetR)

set.seed(42)

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.mergedTranscripts.gtf.tmap"))
tpm_mean_h1299 <- readRDS(file.path(results.dir, "tpm_meanTreatment_from_counts.rds"))

#load aim and cta
ext.dir <- "/omics/groups/OE0219/internal/tinat/external_data/"
AIM <- as.data.frame(data.table::fread(file.path(ext.dir, "Li_etal", "ubiquitous_AIMS.txt")))
CTA <-  as.data.frame(data.table::fread(file.path(ext.dir,  "Cancer_Testis_Antigens.txt"), header=FALSE))

#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c( "DMSO", "SB939","DAC", "DACSB")
comp_col <- col[c(2,8,6)]
names(comp_col)<-c( "DAC_vs_DMSO", "SB939_vs_DMSO",  "DACSB_vs_DMSO")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
ltr_col <- brewer.pal(5, "Accent")
names(ltr_col)<- c("LTR12","LTR12C","LTR12D","LTR12F", "no LTR12")
cell_col  <- brewer.pal(length(unique(colData(vst)$cellline)), "Paired")
names(cell_col)<- as.character(unique(colData(vst)$cellline))

#heatmaps
#define padj cutoff for plotting --> only up
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff
cutoff <- alpha
l2fc <- lfc

#get class code information of de novo assembly, joint with number of exons
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id

#Plot results in Volcano and MA plot {.tabset}
DEG_results_list_plot <- lapply(DEG_results_list, function(x){
  x <- as.data.frame(x)
  #x <- dplyr::left_join(x, anno_classi,by="transcript_id")
  x$col <- "not_significant"
  x$col <- ifelse(test = x$padj < alpha & abs(x$log2FoldChange)>lfc, yes = "significant", 
                  no ="not_significant")
  x$Label <- NA
  x[x$gene_name %in% CTA$V1,][1:5,]$Label <- x[x$gene_name %in% CTA$V1,][1:5,]$gene_name
  x$pvalue <-(-(log10(x$pvalue)))
  x$padj <- (-(log10(x$padj)))
  x
}) 

Volcanos<- list() 
for (i in c( "DAC_vs_DMSO","SB939_vs_DMSO","DACSB_vs_DMSO")){
  Volcanos[[i]] <- ggplot(DEG_results_list_plot[[i]], aes(x=log2FoldChange, y=padj))+
    theme_pubr()+
    geom_point(data = subset(DEG_results_list_plot[[i]], class_code_simple== "non-chimeric (novel)" & col=="not_significant"), 
        fill = class_col[ "non-chimeric (novel)"],color = class_col[ "non-chimeric (novel)"], alpha = 0.2, size=.8) +
    geom_point(data = subset(DEG_results_list_plot[[i]], class_code_simple=="chimeric (novel)" & col=="not_significant"), 
        fill = class_col["chimeric (novel)"], color = class_col["chimeric (novel)"], alpha = 0.2, size=.8) +
    geom_point(data = subset(DEG_results_list_plot[[i]], class_code_simple=="known" & col=="not_significant"), 
        fill = class_col["known"], color = class_col["known"], alpha = 0.2, size=.6) +
    geom_point(data = subset(DEG_results_list_plot[[i]], class_code_simple=="known" & col=="significant"), 
        fill = class_col["known"], color = class_col["known"],alpha = 1, size=.6) +
    geom_point(data = subset(DEG_results_list_plot[[i]], class_code_simple== "chimeric (novel)" & col=="significant"), 
        fill = class_col["chimeric (novel)"],color = class_col["chimeric (novel)"], alpha = 1, size=.8) +
    geom_point(data = subset(DEG_results_list_plot[[i]], class_code_simple== "non-chimeric (novel)" & col=="significant"), 
        fill = class_col[ "non-chimeric (novel)"],color = class_col[ "non-chimeric (novel)"], alpha = 1, size=.8) +
    ggrepel::geom_text_repel(aes(label = Label),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
    geom_vline(xintercept=c(-lfc,lfc), linetype = 2) +
    geom_hline(yintercept=-log10(alpha), linetype = 2)+
    xlab("Transcript expression change\n(log2 fold change)") + ylab("- log10(adj. P value)")
  
  ggsave(plot=Volcanos[[i]],file.path(PostDE.dir,i,"volcano_Class_code_simple_CTALabel.pdf"),height = 5, width = 5, useDingbats = FALSE)
  ggsave(plot=Volcanos[[i]],file.path(PostDE.dir,i,"volcano_Class_code_simple_CTALabel.png"),height = 5, width = 5, device="png")
  print(i)
}


    DEG_results_list_plot[[i]]$class_code_simple <- factor(DEG_results_list_plot[[i]]$class_code_simple, levels=c( "non-chimeric (novel)", "known",  "chimeric (novel)"))
  Volcanos[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="log2FoldChange", y = "padj",
          xlab="Gene expression (log2 fold change)",#xscale="log10" ,
          ylab="- log10(P.adjusted)", #xlim=c(-0.5,0.5),
          color = "class_code_simple" ,fill= "class_code_simple", shape = 16, size = 1.5, palette=class_col, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance",# main=i# Customize reg. line
          label="Label", repel = T, font.label = c( "black"), 
) + rremove("legend") + geom_vline(xintercept=c(-lfc,lfc), linetype = 2) +geom_hline(yintercept=-log10(alpha), linetype = 2)

  pdf(file.path(PostDE.dir,i,"volcano_Class_code_simple_CTALabel.pdf"), height=5, width=5)
  print(Volcanos[[i]])
  dev.off()


 lapply(DEG_results_list[ c( "DAC_vs_DMSO","SB939_vs_DMSO","DACSB_vs_DMSO")], function(x){
    #x <- dplyr::left_join(x, anno_classi,by="transcript_id")
    up <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    print("up")
    print(table(up$class_code_simple))
    down <- x[which(x$padj < cutoff & x$log2FoldChange <(-l2fc)),]
    print("down")
    print(table(down$class_code_simple))
 })