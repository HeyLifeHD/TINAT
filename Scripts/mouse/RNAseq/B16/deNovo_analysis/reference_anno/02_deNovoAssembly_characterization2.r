#start characterizing de novo assembl
#library
library(ggpubr)
library(GenomicRanges)
library(dplyr)
library(ggbreak)
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_analysis"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/visualizations_custom"

#load references
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  import.gff2("/omics/groups/OE0219/internal/genomes/Mmusculus/mm10/rnaseq/gencode.vM19.annotation.gtf")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.mergedTranscripts.gtf.tmap"))
repeats <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/","data","repeats_mm10.rds"))

#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff


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

#plot number of categories$transcript_id <- anno_classi$qry_id
anno_classi_transcript <- unique(anno_classi[, c("transcript_id", "class_code")])
numb_class_code <- as.data.frame(table(anno_classi_transcript$class_code))
colnames(numb_class_code)<- c("class_code", "Freq")
numb_class_code$class_code_simple  <- NA
numb_class_code$class_code_simple <- ifelse(numb_class_code$class_code == "=", "known", NA)
numb_class_code$class_code_simple <- ifelse(numb_class_code$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", numb_class_code$class_code_simple)
numb_class_code$class_code_simple <- ifelse(numb_class_code$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", numb_class_code$class_code_simple)
numb_class_code
numb_class_code$class_code <- factor(numb_class_code$class_code, ordered=TRUE, rev(c("s", "p", "y", "i","x", "u", "c", "o", "m", "n", "k", "j", "=")))
pdf(file.path(output.dir, "transcript_anno_classCode_full_barplot.pdf"))
ggbarplot(numb_class_code,x="Freq", y="class_code", fill="class_code_simple", 
    palette=class_col, xlab="Number of transcripts", ylab="Transcript classification")+ scale_x_break(c(22000, 200000)) 
dev.off()

#visualize transcript classification codes
length(anno_classi$qry_id)
length(unique(anno_classi$qry_id))
table(anno_classi$class_code)
anno_classification <- as.data.frame(table(anno_classi$class_code))
colnames(anno_classification)<- c("class_code", "frequency")
classification_dict_1 <- c("complete,exact match of intron chain","contained in reference (intron compatible)",
    "containment of reference (reverse containment)", "retained intron(s), all introns matched or retained",
    "retained intron(s), not all introns matched/covered", "multi-exon, with at least one junction match",
    "single exon, transfag partially covering an intron, possible pre-mRNA fragment", "other same strand overlap with reference exons",
    "intron match on the opposite strand (likely a mapping error)", "exonic overlap on the opposite strand (like o or e but on the opposite strand",
    "fully contained within a reference intron", " contains a freference within its intron(s)", "possible polymerase run-on (no actual overlap)",
    "repeat (at least 50% bases soft masked)", "none of the above (unknown, intergenic)")
classification_dict_2 <- c("=", "c", "k", "m", "n", "j", "e",
    "o", "s", "x", "i", "y", "p", "r", "u")
classification <- data.frame(class_code=classification_dict_2, explanation=classification_dict_1)
anno_classification <- dplyr::left_join(anno_classification, classification)
p <- ggbarplot(anno_classification,x="class_code", y="frequency", fill="explanation",  label = TRUE, 
    order= anno_classification[order(anno_classification$frequency, decreasing=T),]$class_code,
    ylab="Number of transcripts", legend="right") + 
    rotate()  +
    rremove("ylab")
pdf(file.path(output.dir, "transcript_anno_classCode_barplot.pdf"))
print(p+
    rremove("legend"))
dev.off()
leg <- get_legend(p)
leg <- as_ggplot(leg)
pdf(file.path(output.dir, "transcript_anno_classCode_legend.pdf"))
print(leg)
dev.off()

#annotate repeat anntoated reference as well as DEGs
anno_classi_sub <- anno_classi[, c("class_code","qry_id", "num_exons" )]
colnames(anno_classi_sub)<- c("class_code", "transcript_id", "num_exons")
anno_transcripts <- as.data.frame(anno)[anno$type=="transcript",]
anno_transcripts <- dplyr::left_join(anno_transcripts, anno_classi_sub)
DEG_results_list <- lapply(DEG_results_list, function(x){
    x <- dplyr::left_join(x,anno_classi_sub)
    x$anno_simple <- ifelse(x$class_code == "=", "known", "novel")
    x
})

#export lists
DEG_results_list_exp <- list()
for (i in c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")){
    DEG_results_list_exp[[i]] <- as.data.frame(DEG_results_list[[i]])
    DEG_results_list_exp[[i]]$padj <- round(DEG_results_list_exp[[i]]$padj, 4)
    DEG_results_list_exp[[i]]$log2FoldChange <- round(DEG_results_list_exp[[i]]$log2FoldChange, 4)
    dir.create(file.path(output.dir, i))
    write.table(DEG_results_list_exp[[i]], file=file.path(output.dir,i, "DEG.tsv"),row.names =F, col.names = T , sep="\t", quote=FALSE)
    write.table(DEG_results_list_exp[[i]], file=file.path(output.dir,i, "DEG.csv"),row.names =F, col.names = T , sep=", ", quote=FALSE)
}


#subset upregulated genes
DEG_results_list_up <- lapply(DEG_results_list, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange>lfc),]
    x
})

#look at overlap of novel upregulated
#get transcript ids
DEG_results_list_up_id <- lapply(DEG_results_list_up, function(x){
    x <- x[x$anno_simple=="novel",]
    x$transcript_id
    })
DEG_results_list_up_id <- DEG_results_list_up_id[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
pdf(file.path(output.dir, "Venn_DEG_up_novel.pdf"))
ggVennDiagram::ggVennDiagram(DEG_results_list_up_id, label_alpha = 0, category.names=names(DEG_results_list_up_id))+
    rremove("legend") + 
    theme(text = element_text(size = 8)) #+ 
    #scale_fill_gradient(low="white",high = "#EF8A62")
dev.off()

#get stats of upregulated
DEG_up_stats <- lapply(DEG_results_list_up, function(x){
    x<- x[which(x$padj < 0.01 & x$log2FoldChange>2),]
    x <- as.data.frame(table(x$class_code))
    colnames(x)<- c("class_code", "frequency")
    x <- dplyr::left_join(x, classification)
    x
})
for(i in names(DEG_up_stats)){
    DEG_up_stats[[i]]$comparison <- i
}

#prepare for plotting
DEG_up_stats_sub <- DEG_up_stats[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
DEG_up_stats_sub<- do.call("rbind",DEG_up_stats_sub)

temp <- aggregate(x=DEG_up_stats_sub$frequency, by=list(DEG_up_stats_sub$class_code), FUN=mean)
temp[order(temp$x, decreasing=TRUE),]$Group.1

class_plot <- ggbarplot(DEG_up_stats_sub,x="class_code", y="frequency", fill="comparison",  label = TRUE, 
    ylab="Number of transcripts", order=temp[order(temp$x, decreasing=TRUE),]$Group.1,
    position = position_dodge(0.9), palette=c("#00AFBB", "#E7B800", "#FC4E07")) + 
    rotate() +
    rremove("ylab")
pdf(file.path(output.dir, "DEG_up_anno_classCode_barplot.pdf"))
print(class_plot)
dev.off()

#get number of exons
#get stats of upregulated
DEG_up_stats_exon <- lapply(DEG_results_list_up, function(x){
    x<- x[which(x$padj < 0.01 & x$log2FoldChange>2),]
    x <- x[,c("class_code", "num_exons")]
    x
})
for(i in names(DEG_up_stats_exon)){
    DEG_up_stats_exon[[i]]$comparison <- i
}
#prepare for plotting
DEG_up_stats_exon_sub <- DEG_up_stats_exon[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
DEG_up_stats_exon_sub<- do.call("rbind",DEG_up_stats_exon_sub)
exon_plot <- ggboxplot(DEG_up_stats_exon_sub,x="class_code", y="num_exons", fill="comparison", 
    ylab="Number of exons", order=temp[order(temp$x, decreasing=TRUE),]$Group.1, palette=c("#00AFBB", "#E7B800", "#FC4E07")) + 
    rotate() +
    rremove("ylab")
pdf(file.path(output.dir, "DEG_up_anno_classCode_numExons_boxplot.pdf"))
print(exon_plot)
dev.off()

exon_plot_lim <- ggboxplot(DEG_up_stats_exon_sub,x="class_code", y="num_exons", fill="comparison", outlier.shape=NA,
    ylab="Number of exons", order=temp[order(temp$x, decreasing=TRUE),]$Group.1, palette=c("#00AFBB", "#E7B800", "#FC4E07")) + 
    rremove("ylab")+ coord_flip()+ylim(c(0,40))
pdf(file.path(output.dir, "DEG_up_anno_classCode_numExons_boxplot_limited.pdf"))
print(exon_plot_lim)
dev.off()

#combined plot
pdf(file.path(output.dir, "DEG_up_anno_classCode_combined.pdf"))
ggarrange(class_plot, exon_plot_lim, common.legend = TRUE,
          ncol = 2, nrow = 1)
dev.off()

#plot distance to ltr 12 of annotated and non annotated
for(i in names(DEG_results_list_up)){
    DEG_results_list_up[[i]]$comparison <- i
}
DEG_results_list_up_sub <- DEG_results_list_up[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
DEG_results_list_up_sub_unlist <- do.call("rbind",DEG_results_list_up_sub)
pdf(file.path(output.dir, "DEG_up_anno_distLTR12_histo.pdf"))
gghistogram(DEG_results_list_up_sub_unlist, x="dist_nearest_LTR12repeat", 
    fill="anno_simple", facet.by="comparison", add="mean")
dev.off()


#plot known and novel transcripts of upregualted degs
#get stats of upregulated
DEG_up_stats <- lapply(DEG_results_list_up, function(x){
    x<- x[which(x$padj < 0.01 & x$log2FoldChange>2),]
    x <- as.data.frame(table(x$anno_simple))
    colnames(x)<- c("transcript_class", "frequency")
    x
})
#prepare for plotting
DEG_up_stats_sub <- DEG_up_stats[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
#DEG_up_stats_sub<- do.call("rbind",DEG_up_stats_sub)
for(i in names(DEG_up_stats_sub)){
    labs <- paste0(DEG_up_stats_sub[[i]]$frequency, " ",DEG_up_stats_sub[[i]]$transcript_class," transcripts")
    dir.create(file.path(output.dir, i))
    pdf(file.path(output.dir, i, "DEG_up_anno_simple.pdf"))
    print(ggpie(DEG_up_stats_sub[[i]], "frequency", label = labs,
        lab.pos = "in", lab.font = "white",
        fill = "transcript_class", color = "white"))
    dev.off()
}

#annotate rep class of dEG up
for(i in names( DEG_results_list)){
    DEG_results_list[[i]]$ERE_TSS_anno <- ifelse(DEG_results_list[[i]]$dist_nearest_repeat == 0, DEG_results_list[[i]]$nearest_repeat_repClass, "noERE")
    DEG_results_list[[i]]$ERE_TSS_anno <- ifelse(DEG_results_list[[i]]$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "noERE"), DEG_results_list[[i]]$ERE_TSS_anno , "other")
}
#get stats of upregulated
DEG_up_stats <- lapply(DEG_results_list, function(x){
    x<- x[which(x$padj < 0.01 & x$log2FoldChange>2),]
    x <- as.data.frame(table(x$ERE_TSS_anno))
    colnames(x)<- c("repeatClass", "frequency")
    x
})
#prepare for plotting
DEG_up_stats_sub <- DEG_up_stats[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
#DEG_up_stats_sub<- do.call("rbind",DEG_up_stats_sub)
for(i in names(DEG_up_stats_sub)){
    labs <- paste0(DEG_up_stats_sub[[i]]$frequency, " ",DEG_up_stats_sub[[i]]$repeatClass," transcripts")
    pdf(file.path(output.dir, i, "DEG_up_anno_repClass.pdf"))
    print(ggpie(DEG_up_stats_sub[[i]], "frequency", label = labs,
        lab.pos = "in", lab.font = "white",
        fill = "repeatClass", color = "white"))
    dev.off()
}

#same stratified in novel and known
DEG_results_list_up <- list()
DEG_results_list_up_split <-list()
for(i in  c("DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")){
    DEG_results_list_up[[i]]<- DEG_results_list[[i]][which(DEG_results_list[[i]]$padj < 0.01 & DEG_results_list[[i]]$log2FoldChange>2),]
    DEG_results_list_up_split[[i]]<- split(DEG_results_list_up[[i]], DEG_results_list_up[[i]]$anno_simple)
    DEG_up_stats_sub <- lapply(DEG_results_list_up_split[[i]], function(x){
            x <- as.data.frame(table(x$ERE_TSS_anno))
            colnames(x)<- c("repeatClass", "frequency")
            x
    })
    labs <- paste0(DEG_up_stats_sub[["known"]]$frequency, " ",DEG_up_stats_sub[["known"]]$repeatClass," transcripts")
    known <- ggpie(DEG_up_stats_sub[["known"]], "frequency", label = labs,
        lab.pos = "in", lab.font = "white",
        fill = "repeatClass", color = "white")
    labs <- paste0(DEG_up_stats_sub[["novel"]]$frequency, " ",DEG_up_stats_sub[["novel"]]$repeatClass," transcripts")
    novel <- ggpie(DEG_up_stats_sub[["novel"]], "frequency", label = labs,
        lab.pos = "in", lab.font = "white",
        fill = "repeatClass", color = "white")
    pdf(file.path(output.dir, i, "DEG_up_anno_repClass_strat_simpleanno.pdf"))
        print(ggarrange(known, novel,  common.legend = TRUE,
          labels = c("known", "novel"),
          ncol = 2, nrow = 1
          ))
    dev.off()
}