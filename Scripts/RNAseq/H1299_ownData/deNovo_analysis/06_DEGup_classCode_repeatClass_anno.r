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
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
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
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))

#Take a look at design and annotation of samples
design(dds)
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff


#heatmaps
#define padj cutoff for plotting --> only up
cutoff <- alpha
l2fc <- lfc
DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
genes2plot <- lapply(DEG_results_list_sub, function(x){
  x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
  x <- rownames(x)
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]
lapply(genes2plot, function(x)length(x))

#heatmap
mat_genes<- assay(vst)
col <- c("#00AFBB", "#E7B800", "#FC4E07","gray")
names(col) <- levels(colData(vst)$treatment)
anno_colors <- list(treatment=col)
anno <- colData(vst)
annovst <- anno[,c("treatment" ), drop=FALSE]
DEG <- unique(unlist(genes2plot))

length(DEG)
plot <- mat_genes[which(rownames(mat_genes) %in% DEG ),]
p <- pheatmap(plot, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),fontsize=16,
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
                       #cutree_rows =5, 
                       file=file.path(PostDE.dir, "Heat_upDEG_subComparison_correlation.pdf")
                       ) 


#look at overlap of upregulated degs and transcript class annotation
#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi_split <- split(anno_classi, anno_classi$class_code_simple)
anno_classi_split <- lapply(anno_classi_split, function(x){unique(x$qry_id)})
compar <- c(genes2plot,anno_classi_split )
#plot
col <- c("#00AFBB", "#E7B800", "#FC4E07", "orange", "tomato3")#"darkgray"
names(col)<- c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO", "chimeric (novel)", "non-chimeric (novel)")# "known", 
#col <- c( "#FC4E07","#00AFBB", "tomato3", "#FC4E07", "orange","darkgray")

pdf(file.path(PostDE.dir, "Upset_upDEG_transcript_class_simple.pdf"), width=7.5, height=5)
UpSetR::upset(UpSetR::fromList(compar[-5]),  order.by = "freq", set_size.show = TRUE ,nsets=6, sets.bar.color=col)
dev.off()
###to do --> modify source code to not scale data and only modify axis

#plot class code simple of transcripts of upregualted degs with erg distance 
#get transcript and repeat annotation
anno_classi$transcript_id <- anno_classi$qry_id 
DEG_up_anno <- lapply(DEG_results_list_sub, function(x){
    x<- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x$ERE_TSS_anno <- NA
    x$ERE_TSS_anno <- ifelse(x$dist_nearest_repeat == 0, x$nearest_repeat_repClass, "no ERE")
    x$ERE_TSS_anno <- ifelse(x$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), x$ERE_TSS_anno , "other")
    x
})
#split by anno simple and get stats
DEG_up_anno_split_stats <- lapply(DEG_up_anno, function(x){
    x <- split(x, x$class_code_simple)
    x <- lapply(x, function(y){
    y <- as.data.frame(table(y$ERE_TSS_anno))
    colnames(y) <- c("repeat_class", "Freq")
    y
    })
    #x <- do.call("rbind",x)
    #x$class_code_simple <- sapply(strsplit(rownames(x), ".", fixed=TRUE), "[",1)
    #x
})
#plot
pal <- c("deeppink3", "red4", "midnightblue", "gold1", "salmon4")
names(pal)<-  levels(DEG_up_anno_split_stats$SB939_vs_DMSO$repeat_class)  

pie <- list()
for (i in names(DEG_up_anno_split_stats)){
    pie[[i]]<- list()
    for(j in names(DEG_up_anno_split_stats[[i]])){
        labs <- paste0(DEG_up_anno_split_stats[[i]][[j]]$Freq)
        pie[[i]][[j]] <- ggpie(DEG_up_anno_split_stats[[i]][[j]], "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=pal,
                fill = "repeat_class", color = "black",title=j)+rremove("legend.title")
    }
}

pdf(file.path(PostDE.dir, "Pie_upDEG_transcript_class_ERE_anno.pdf"), width=10, height=10)
    ggarrange(pie$"DAC_vs_DMSO"$"non-chimeric (novel)",pie$"SB939_vs_DMSO"$"non-chimeric (novel)", pie$"DACandSB939_vs_DMSO"$"non-chimeric (novel)", 
        pie$"DAC_vs_DMSO"$"chimeric (novel)",pie$"SB939_vs_DMSO"$"chimeric (novel)", pie$"DACandSB939_vs_DMSO"$"chimeric (novel)",  
        pie$"DAC_vs_DMSO"$"known",pie$"SB939_vs_DMSO"$"known", pie$"DACandSB939_vs_DMSO"$"known",   
        common.legend = TRUE,
        labels = names(pie),
        ncol = 3, nrow = 3
        )
dev.off()
pdf(file.path(PostDE.dir, "Pie_upDEG_transcript_class_ERE_anno_differentOrder.pdf"), width=10, height=10)
    ggarrange(pie$"DACandSB939_vs_DMSO"$"non-chimeric (novel)",  pie$"DACandSB939_vs_DMSO"$"chimeric (novel)", pie$"DACandSB939_vs_DMSO"$"known",  
        pie$"DAC_vs_DMSO"$"non-chimeric (novel)", pie$"DAC_vs_DMSO"$"chimeric (novel)", pie$"DAC_vs_DMSO"$"known",
        pie$"SB939_vs_DMSO"$"non-chimeric (novel)",    pie$"SB939_vs_DMSO"$"chimeric (novel)", pie$"SB939_vs_DMSO"$"known",   
        common.legend = TRUE,
        labels = names(pie[[1]]),
        ncol = 3, nrow = 3
        )
dev.off()


pdf(file.path(PostDE.dir, "Pie_upDEG_transcript_class_ERE_anno_differentOrder_sub.pdf"), width=5, height=10)
    ggarrange(pie$"SB939_vs_DMSO"$"non-chimeric (novel)",    pie$"SB939_vs_DMSO"$"chimeric (novel)", 
        pie$"DAC_vs_DMSO"$"non-chimeric (novel)", pie$"DAC_vs_DMSO"$"chimeric (novel)",
        pie$"DACandSB939_vs_DMSO"$"non-chimeric (novel)",  pie$"DACandSB939_vs_DMSO"$"chimeric (novel)", 
        common.legend = TRUE,
        labels = c("non-chimeric (novel)", "chimeric (novel)"),
        ncol = 2, nrow = 3
        )
dev.off()


#focus on novel upregulated transcripts
anno_classi$transcript_id <- anno_classi$qry_id
DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
DEG_results_list_sub_anno <- lapply(DEG_results_list_sub, function(x){
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x
})
DEG_results_list_sub_anno_sub <- lapply(DEG_results_list_sub_anno, function(x){
      x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
      x <- x[x$class_code_simple != "known",]
      print(nrow(x))
      x
}) 

#plot as an heatmap
genes2plot <- lapply(DEG_results_list_sub_anno_sub, function(x){
  x <-x$transcript_id
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]
lapply(genes2plot, function(x)length(x))

#heatmap
mat_genes<- assay(vst)
DEG <- unique(unlist(genes2plot))
length(DEG)
plot <- mat_genes[which(rownames(mat_genes) %in% DEG ),]
p <- pheatmap(plot, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),fontsize=16,
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
                       #cutree_rows =5, 
                       file=file.path(PostDE.dir, "Heat_upDEG_novel_subComparison_correlation.pdf")
                       ) 
p <- pheatmap(t(plot), scale="column", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_row=as.data.frame(annovst),fontsize=16,
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
                       #cutree_rows =5, 
                       file=file.path(PostDE.dir, "Heat_upDEG_novel_subComparison_correlation_t.pdf")
                       ) 
#changed order
p <- pheatmap(plot[,c(1,5,9,3,7,11,2,6,10, 4,8,12)], scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),fontsize=16,cluster_cols=FALSE,
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
                       file=file.path(PostDE.dir, "Heat_upDEG_novel_subComparison_correlation_noColCluster.pdf")
                       ) 
p <- pheatmap(t(plot[,c(1,5,9,3,7,11,2,6,10, 4,8,12)]), scale="column", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_row=as.data.frame(annovst),fontsize=16,cluster_rows=FALSE,
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_cols="correlation" ,
                       file=file.path(PostDE.dir, "Heat_upDEG_novel_subComparison_correlation_noColCluster_T.pdf"), width=7
                       ) 
pdf(file.path(PostDE.dir, "Venn_upDEG_transcript_class_code_simple_novel.pdf"), height=5, width=5)
ggVennDiagram::ggVennDiagram(genes2plot[c(2,1,3)]) + scale_fill_gradient(low="white",high =  rev(hcl.colors(10,"Reds"))[9])
dev.off()

#venn diagramm of upregulated degs with novel genes
DEG_results_list_sub_anno_sub <- lapply(DEG_results_list_sub_anno, function(x){
      x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
      #x <- x[x$class_code_simple != "known",]
      print(nrow(x))
      x
}) 
 
#plot as an heatmap
genes2plot <- lapply(DEG_results_list_sub_anno_sub, function(x){
  x <-x$transcript_id
  x
})
anno_classi_split_sub <- anno_classi_split[-2]
genes2plot <- c(genes2plot,list(novel=unlist(anno_classi_split_sub)))
pdf(file.path(PostDE.dir, "Venn_upDEG_transcript_class_code_simple.pdf"))
ggVennDiagram::ggVennDiagram(genes2plot) + scale_fill_gradient(low="white",high = hcl.colors(10,"Reds")[10])
dev.off()

#focus on combined novel class
#plot class code simple of transcripts of upregualted degs with erg distance 
#get transcript and repeat annotation
anno_classi$transcript_id <- anno_classi$qry_id 
DEG_up_anno <- lapply(DEG_results_list_sub, function(x){
    x<- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x$ERE_TSS_anno <- NA
    x$ERE_TSS_anno <- ifelse(x$dist_nearest_repeat == 0, x$nearest_repeat_repClass, "no ERE")
    x$ERE_TSS_anno <- ifelse(x$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), x$ERE_TSS_anno , "other")
    x
})
#split by anno simple and get stats
DEG_up_anno_split_stats <- lapply(DEG_up_anno, function(x){
    x <- x[x$class_code_simple != "known",]
    x <- as.data.frame(table(x$ERE_TSS_anno))
    colnames(x) <- c("repeat_class", "Freq")
    x

})
#plot
pal <- c("deeppink3", "red4", "midnightblue", "gold1", "salmon4")
names(pal)<-  levels(DEG_up_anno_split_stats$SB939_vs_DMSO$repeat_class)  

pie <- list()
for (i in names(DEG_up_anno_split_stats)){
        labs <- paste0(DEG_up_anno_split_stats[[i]]$Freq)
        pie[[i]] <- ggpie(DEG_up_anno_split_stats[[i]], "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=pal,
                fill = "repeat_class", color = "black",title=j)+rremove("legend.title")
 
}

pdf(file.path(PostDE.dir, "Pie_upDEG_transcript_class_ERE_anno_onlyNovel.pdf"), width=5, height=10)
    ggarrange(pie$"SB939_vs_DMSO", pie$"DAC_vs_DMSO",pie$"DACandSB939_vs_DMSO",
        common.legend = TRUE,
        labels = names(pie),
        ncol = 1, nrow = 3
        )
dev.off()