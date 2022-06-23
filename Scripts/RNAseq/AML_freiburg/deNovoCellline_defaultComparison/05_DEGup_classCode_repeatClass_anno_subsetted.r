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
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220622_AMLPatient_deNovoCellline_quantification_analysis_revised/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE_subsetted")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

anno <- readRDS("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.mergedTranscripts.gtf.tmap"))

#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2)]
#names(treat_col)<-c( "DMSO", "SB939","DAC", "DACSB")
names(treat_col)<-c( "scr", "c1d1","c1d8", "c1d9")
#comp_col <- col[c(2,8,6)]
#names(comp_col)<-c( "DAC_vs_DMSO", "SB939_vs_DMSO",  "DACandSB939_vs_DMSO")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
ltr_col <- brewer.pal(5, "Accent")
names(ltr_col)<- c("LTR12","LTR12C","LTR12D","LTR12F", "no LTR12")
#cell_col  <- brewer.pal(length(unique(colData(vst)$Patient_ID)), "Paired")
#names(cell_col)<- as.character(unique(colData(vst)$Patient_ID))
cell_col  <- c(brewer.pal(10, "Paired"), brewer.pal(4, "Dark2"))
names(cell_col)<- as.character(unique(colData(vst)$Patient_ID))
font_size=12

#heatmaps
#define padj cutoff for plotting --> only up
alpha <- 0.05 #set FDR cutoff
lfc <- 1##set logfold2 cutoff

cutoff <- alpha
l2fc <- lfc
DEG_results_list_sub <- DEG_results_list
genes2plot <- lapply(DEG_results_list_sub, function(x){
  x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
  x <- rownames(x)
  print(length(x))
  x
})

# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]
lapply(genes2plot, function(x)length(x))

genes2plot$combinedDEGs <-  unique(unlist(genes2plot))
#get class code information of de novo assembly, joint with number of exons
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id

for( i in names(genes2plot)){
    dir.create(file.path(PostDE.dir,i))
    #subsetgenes of interest
    genes_oi <- genes2plot[[i]]
    anno_classi_sub <- anno_classi[anno_classi$transcript_id %in% genes_oi,]
    anno_sub <- anno[anno$transcript_id%in% genes_oi, ]
    #rename class codes
    print(table(anno_classi_sub$class_code_simple))
    trans_class <- as.data.frame(table(anno_classi_sub$class_code_simple))

    #get transcript numbers
    class_plot <- ggbarplot(trans_class, x="Var1", y="Freq",
        fill="Var1",  palette=class_col,yscale="log10",lab.pos="out",#color="Var1",
        ylab="Number of transcripts", label=TRUE, order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
        font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
        )    +
        ggpubr::rotate()+
        rremove("legend") +
        rremove("ylab")
    pdf(file.path(PostDE.dir,i,  "transcript_classification_simplified.pdf"), height=5, width=5 )
    print(class_plot)
    dev.off()

    labs <- paste0(trans_class$Freq)
    class_plot_pie <- ggpie(trans_class, "Freq", label = labs,legend="right",
            lab.pos = "in", lab.font = "black",label.size = 16, palette=class_col,
            fill = "Var1", color = "black")+rremove("legend.title")
    pdf(file.path(PostDE.dir,i, "transcript_classification_simplified_pie.pdf"), height=5, width=5)
    print(class_plot_pie)
    dev.off()

    exon_plot <- ggboxplot(anno_classi_sub, x="class_code_simple", y="num_exons",
        fill="class_code_simple",  palette=class_col,#color="class_code_simple",
        ylab="Number of exons", order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
        outlier.shape=NA,
        font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
        )    +
        rremove("legend") +
        rremove("ylab") +
        coord_flip() +
        ylim(c(0,40))
    pdf(file.path(PostDE.dir,i, "transcript_exon_number_simplified.pdf"), height=5, width=5)
    print(exon_plot)
    dev.off()
    

    #plot mono or multi exonic in pie chart
    anno_classi_sub$exon_class <- ifelse(anno_classi_sub$num_exons==1, "mono-exonic", "multi-exonic")
    anno_classi_split <- split(anno_classi_sub, anno_classi_sub$class_code_simple)
    exon_stats <- lapply(anno_classi_split, function(x){
            x <-  as.data.frame(table(x$exon_class))

    })
    exon_pie <- list()
    for(m in names(exon_stats)){
         labs <- paste0(exon_stats[[m]]$Freq)
        exon_pie[[m]]<- ggpie(exon_stats[[m]], "Freq", label = labs,
                    lab.pos = "in", lab.font = "black", palette=exon_col,#label.size = 16,
                fill = "Var1", color = "black")+rremove("legend.title")
    }
    pies <- ggarrange(exon_pie[[trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1[3]]],
                exon_pie[[trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1[2]]], exon_pie[[trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1[1]]],
                 common.legend = TRUE,ncol = 1, nrow = 3
            )
    pdf(file.path(PostDE.dir,i, "Exon_class.pdf"), width=5)
        print(pies)
    dev.off()

    #get length of transcripts
    transcripts_width <- data.frame(transcript_id=anno_sub[anno_sub$type=="transcript",]$transcript_id, width=width(anno_sub[anno_sub$type=="transcript",]))
    anno_classi_sub$transcript_id <- anno_classi_sub$qry_id
    anno_classi_sub <- dplyr::left_join(anno_classi_sub, transcripts_width, by="transcript_id")
    anno_classi_sub$width <- anno_classi_sub$width/1000

    length_plot <- ggboxplot(anno_classi_sub, x="class_code_simple", y="width",
        fill="class_code_simple",  palette=class_col,#color="class_code_simple",
        ylab="Length of transcripts [kbp]", order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
        outlier.shape=NA,
        font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
        )    +
        rremove("legend") +
        rremove("ylab") +
        coord_flip() +#
    ylim(c(0,200))
    pdf(file.path(PostDE.dir,i, "transcript_length_simplified.pdf"), height=5, width=5)
    print(length_plot)
    dev.off()

    #combined plot
    pdf(file.path(PostDE.dir,i, "transcript_anno_classCode_simple_length_exonNum_combined.pdf"), width=10, height=5)
    print(ggarrange(class_plot, length_plot+rremove("y.text"), exon_plot+rremove("y.text"), common.legend = FALSE,
            ncol = 3, nrow = 1, widths=c(1.5,1,1)
            ))
    dev.off()
    pdf(file.path(PostDE.dir,i, "transcript_anno_classCodePie_simple_length_exonNum_combined.pdf"), width=10, height=5)
    print(ggarrange(class_plot_pie+rremove("legend"), length_plot, exon_plot+rremove("y.text"), common.legend = FALSE,
            ncol = 3, nrow = 1, widths=c(1,1.5,1)
            ))
    dev.off()


    #plot class code simple of transcripts of upregualted degs with erg distance 
    #get transcript and repeat annotation
    anno_classi_sub$transcript_id <- anno_classi_sub$qry_id 
    anno_sub_transcript <- as.data.frame(anno_sub[anno_sub$type=="transcript",])
    anno_classi_ere <- dplyr::left_join(anno_classi_sub, anno_sub_transcript, by="transcript_id")
    anno_classi_ere$ERE_TSS_anno <- NA
    anno_classi_ere$ERE_TSS_anno <- ifelse(anno_classi_ere$dist_nearest_repeat == 0, anno_classi_ere$nearest_repeat_repClass, "no ERE")
    anno_classi_ere$ERE_TSS_anno <- ifelse(anno_classi_ere$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), anno_classi_ere$ERE_TSS_anno , "other")

    #split
    anno_classi_ere_split <- split(anno_classi_ere, anno_classi_ere$class_code_simple)
    anno_classi_ere_split_stats <- lapply(anno_classi_ere_split, function(y){
    y <- as.data.frame(table(y$ERE_TSS_anno))
    colnames(y) <- c("repeat_class", "Freq")
    y
    })

    #plot
    pie_2 <- list()
        for(j in names(anno_classi_ere_split_stats)){
            labs <- paste0(anno_classi_ere_split_stats[[j]]$Freq)
            pie_2[[j]] <- ggpie(anno_classi_ere_split_stats[[j]], "Freq", #label = labs,
                    lab.pos = "in", lab.font = "black", palette=ere_col,#label.size = 16,
                    fill = "repeat_class", color = "black")+rremove("legend.title")
        }

    pie_eres <- ggarrange(pie_2[[trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1[3]]],
        pie_2[[trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1[2]]], pie_2[[trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1[1]]],
        labels=NULL,
        common.legend = TRUE,
        ncol = 1, nrow = 3
        )
    pdf(file.path(PostDE.dir,i, "Pie_transcript_class_ERE_anno.pdf"), width=10, height=10)
    print(pie_eres)
    dev.off()


    pdf(file.path(PostDE.dir, i,"transcript_anno_classCodePie_simple_length_exonNum_exonPies_ere_anno_combined_v2.pdf"), width=15, height=5)
    print(ggarrange(class_plot_pie+rremove("legend"), length_plot, exon_plot+rremove("y.text"),pies,pie_eres, common.legend = FALSE,
            ncol = 5, nrow = 1, widths=c(1,2.5,1.5,1,1),labels=c("a", "b", "c", "d","e"))
            )
    dev.off() 
    
    # #plot tpm
    # tpm_mean_h1299 <- as.data.frame(tpm_mean_h1299)
    # tpm_mean_h1299$transcript_id <- rownames(tpm_mean_h1299)
    # tpm_mean_h1299_sub <- tpm_mean_h1299[tpm_mean_h1299$transcript_id %in% genes_oi,]
    # tpm_mean_h1299_merged <- dplyr::left_join(tpm_mean_h1299_sub, anno_classi[, c("transcript_id", "class_code_simple")])
    # tpm_mean_h1299_merged$transcript_id <- NULL
    # tpm_mean_h1299_merged_melt <- reshape2::melt(tpm_mean_h1299_merged)
    # tpm_mean_h1299_merged_melt$class_code_simple<- factor(tpm_mean_h1299_merged_melt$class_code_simple, 
    #     levels=   trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1)

     
    # pal <- c("darkgray", "orange", "tomato3")
    # font_size <- c(16,"plain", "black")
    # names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
    # tpm_allTranscripts <- ggboxplot(tpm_mean_h1299_merged_melt, x="variable", y="value",
    #     fill="class_code_simple",  palette=class_col,lab.pos="out",#color="Var1",
    #     ylab="TPM",  order=rev(c("DMSO", "SB939", "DAC", "DACSB")),outlier.shape=NA,
    #     font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size, 
    #     ) +
    #     rremove("legend.title") +
    #     rremove("ylab") +
    #     coord_flip() +and
    #     ylim(c(0,10))
    # pdf(file.path(PostDE.dir, i,"transcript_classification_TPM.pdf"), height=5, width=5 )
    # print(tpm_allTranscripts)
    # dev.off()
}