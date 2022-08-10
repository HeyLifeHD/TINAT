#start characterizing de novo assembly compared to reference
#library
library(ggpubr)
library(GenomicRanges)
library(dplyr)

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/visualizations_custom"

#load references
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))

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
glob_col <-  c("gray", col[c(10)])
names(glob_col)<- c("genecode", "deNovo")

#get numbers of novel genes and transcripts
numb_trans_deNovo <- length(unique(anno_classi$qry_id))
numb_genes_deNovo <- length(unique(anno_classi$qry_gene_id))
numb_trans_genecode <- length(unique(anno_original$transcript_id))
numb_genes_genecode<- length(unique(anno_original$gene_name))
gene_trans_numbers <- data.frame(assembly=c("genecode","deNovo"),
    genes=c(numb_genes_genecode, numb_genes_deNovo),
    transcripts=c(numb_trans_genecode, numb_trans_deNovo))
gene_trans_numbers <- reshape2::melt(gene_trans_numbers)
gene_trans_numbers$assembly <- factor(gene_trans_numbers$assembly, ordered = TRUE, 
levels = c("genecode", "deNovo"))

#plotting
font_size <- c(16,"plain", "black")
pdf(file.path(output.dir, "number_genes_transcripts_log10.pdf"), height=5, width=5)
ggbarplot(gene_trans_numbers, x="variable", y="value",
    fill="assembly", palette=glob_col,yscale="log10",#color="assembly", 
    ylab="Number", label=TRUE, position = position_dodge(0.8),
    legend.title="Assembly",
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    )    +
    rremove("legend.title") +
    rremove("xlab")
dev.off()
pdf(file.path(output.dir, "number_genes_transcripts.pdf"), height=5, width=5)
ggbarplot(gene_trans_numbers, x="variable", y="value",
    fill="assembly", palette=glob_col,# color="assembly",
    ylab="Number", label=TRUE, position = position_dodge(0.8),
    legend.title="Assembly",
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    )    +
    rremove("legend.title") +
    rremove("xlab")
dev.off()



#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
table(anno_classi$class_code_simple)
trans_class <- as.data.frame(table(anno_classi$class_code_simple))

#get transcript numbers
class_plot <- ggbarplot(trans_class, x="Var1", y="Freq",
    fill="Var1",  palette=class_col,yscale="log10",lab.pos="out",#color="Var1",
    ylab="Number of transcripts", label=TRUE, order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    )    +
    ggpubr::rotate()+
    rremove("legend") +
    rremove("ylab")
pdf(file.path(output.dir, "transcript_classification_simplified.pdf"), height=5, width=5 )
print(class_plot)
dev.off()

labs <- paste0(trans_class$Freq)
class_plot_pie <- ggpie(trans_class, "Freq", label = labs,legend="right",
        lab.pos = "in", lab.font = "black",label.size = 16, palette=class_col,
        fill = "Var1", color = "black")+rremove("legend.title")
pdf(file.path(output.dir, "transcript_classification_simplified_pie.pdf"), height=5, width=5)
print(class_plot_pie)
dev.off()

#get number of exons
exon_plot <- ggboxplot(anno_classi, x="class_code_simple", y="num_exons",
    fill="class_code_simple",  palette=class_col,#color="class_code_simple",
    ylab="Number of exons", order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
    outlier.shape=NA,
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    )    +
    rremove("legend") +
    rremove("ylab") +
    coord_flip() +
    ylim(c(0,40))
pdf(file.path(output.dir, "transcript_exon_number_simplified.pdf"), height=5, width=5)
print(exon_plot)
dev.off()


#get length of transcripts
transcripts_width <- data.frame(transcript_id=anno[anno$type=="transcript",]$transcript_id, width=width(anno[anno$type=="transcript",]))
anno_classi$transcript_id <- anno_classi$qry_id
anno_classi <- dplyr::left_join(anno_classi, transcripts_width, by="transcript_id")
anno_classi$width <- anno_classi$width/1000

length_plot <- ggboxplot(anno_classi, x="class_code_simple", y="width",
    fill="class_code_simple",  palette=class_col,#color="class_code_simple",
    ylab="Length of transcripts [kbp]", order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
    outlier.shape=NA,
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    )    +
    rremove("legend") +
    rremove("ylab") +
    coord_flip() +#
   ylim(c(0,200))
pdf(file.path(output.dir, "transcript_length_simplified.pdf"), height=5, width=5)
print(length_plot)
dev.off()

#combined plot
pdf(file.path(output.dir, "transcript_anno_classCode_simple_length_exonNum_combined.pdf"), width=10, height=5)
ggarrange(class_plot, length_plot+rremove("y.text"), exon_plot+rremove("y.text"), common.legend = FALSE,
          ncol = 3, nrow = 1, widths=c(1.5,1,1)
          )
dev.off()
pdf(file.path(output.dir, "transcript_anno_classCodePie_simple_length_exonNum_combined.pdf"), width=10, height=5)
ggarrange(class_plot_pie+rremove("legend"), length_plot, exon_plot+rremove("y.text"), common.legend = FALSE,
          ncol = 3, nrow = 1, widths=c(1,1.5,1)
          )
dev.off()

#plot mono or multi exonic in pie chart
anno_classi$exon_class <- ifelse(anno_classi$num_exons==1, "mono-exonic", "multi-exonic")
anno_classi_split <- split(anno_classi, anno_classi$class_code_simple)
exon_class <- lapply(anno_classi_split, function(x){
    as.data.frame(table(x$exon_class))
})

labs <- paste0(exon_class[["non-chimeric (novel)"]]$Freq)
non_chimeric <- ggpie(exon_class[["non-chimeric (novel)"]], "Freq", label = labs,
        lab.pos = "in", lab.font = "black", palette=exon_col,#label.size = 16,
        fill = "Var1", color = "black")+rremove("legend.title")
labs <- paste0(exon_class[["chimeric (novel)"]]$Freq)
chimeric <- ggpie(exon_class[["chimeric (novel)"]], "Freq", label = labs,
        lab.pos = "in", lab.font = "black", palette=exon_col,#label.size = 16,
        fill = "Var1", color = "black")+rremove("legend.title")
labs <- paste0(exon_class[["known"]]$Freq)
known <- ggpie(exon_class[["known"]], "Freq", label = labs,
        lab.pos = "in", lab.font = "black", palette=exon_col,#label.size = 16,
        fill = "Var1", color = "black")+rremove("legend.title")
pies <- ggarrange(non_chimeric, chimeric,known,  common.legend = TRUE,
          #labels = c("non-chimeric (novel)", "chimeric (novel)", "known"),
          ncol = 1, nrow = 3
          )
pdf(file.path(output.dir, "Exon_class.pdf"), width=5)
    print(pies)
dev.off()

pdf(file.path(output.dir, "transcript_anno_classCodePie_simple_length_exonNum_exonPies_combined.pdf"), width=15, height=5)
ggarrange(class_plot_pie+rremove("legend"), length_plot, exon_plot+rremove("y.text"),pies, common.legend = FALSE,
          ncol = 4, nrow = 1, widths=c(1,2.5,1.5,1)
          )
dev.off()

pdf(file.path(output.dir, "transcript_anno_classCodePie_simple_length_exonNum_exonPies_combined_v2.pdf"), width=15, height=5)
ggarrange(class_plot_pie+rremove("legend"), length_plot, exon_plot+rremove("y.text"),pies, common.legend = FALSE,
          ncol = 4, nrow = 1, widths=c(1,2.5,1.5,1),labels=c("a", "b", "c", "d")
          )
dev.off() 


#plot class code simple of transcripts of upregualted degs with erg distance 
#get transcript and repeat annotation
anno_classi$transcript_id <- anno_classi$qry_id 
anno_transcript <- as.data.frame(anno[anno$type=="transcript",])
anno_classi_ere <- dplyr::left_join(anno_classi, anno_transcript, by="transcript_id")
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
pal <- c("deeppink3", "red4", "midnightblue", "gold1", "salmon4")
names(pal)<-  levels(anno_classi_ere_split_stats$repeat_class)  

pie_ere <- list()
    for(j in names(anno_classi_ere_split_stats)){
        labs <- paste0(anno_classi_ere_split_stats[[j]]$Freq)
        pie[[j]] <- ggpie(anno_classi_ere_split_stats[[j]], "Freq", #label = labs,
                lab.pos = "in", lab.font = "black", palette=ere_col,#label.size = 16,
                fill = "repeat_class", color = "black")+rremove("legend.title")
    }

pie_eres <- ggarrange(pie$"non-chimeric (novel)",
        pie$"chimeric (novel)",
        pie$"known",labels=NULL,
        common.legend = TRUE,
        ncol = 1, nrow = 3
        )
pdf(file.path(output.dir, "Pie_transcript_class_ERE_anno.pdf"), width=10, height=10)
pie_eres
dev.off()


pdf(file.path(output.dir, "transcript_anno_classCodePie_simple_length_exonNum_exonPies_ere_anno_combined_v2.pdf"), width=15, height=5)
ggarrange(class_plot_pie+rremove("legend"), length_plot, exon_plot+rremove("y.text"),pies,pie_eres, common.legend = FALSE,
          ncol = 5, nrow = 1, widths=c(1,2.5,1.5,1,1),labels=c("a", "b", "c", "d","e")
          )
dev.off() 