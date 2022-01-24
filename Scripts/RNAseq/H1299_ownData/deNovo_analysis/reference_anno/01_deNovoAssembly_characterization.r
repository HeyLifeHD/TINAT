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
pal <- c("darkgray", "orange")
names(pal)<- c("genecode", "deNovo")
font_size <- c(16,"plain", "black")
pdf(file.path(output.dir, "number_genes_transcripts.pdf"), height=7, width=7)
ggbarplot(gene_trans_numbers, x="variable", y="value",
    fill="assembly", color="assembly", palette=pal,yscale="log10",
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
pal <- c("darkgray", "orange", "tomato3")
names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
class_plot <- ggbarplot(trans_class, x="Var1", y="Freq",
    fill="Var1", color="Var1", palette=pal,yscale="log10",lab.pos="out",
    ylab="Number of transcripts", label=TRUE, order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    )    +
    rotate()+
    rremove("legend") +
    rremove("ylab")
pdf(file.path(output.dir, "transcript_classification_simplified.pdf"))
print(class_plot)
dev.off()

#get number of exons
exon_plot <- ggboxplot(anno_classi, x="class_code_simple", y="num_exons",
    fill="class_code_simple", color="class_code_simple", palette=pal,
    ylab="Number of exons", order=trans_class[order(trans_class$Freq, decreasing=TRUE),]$Var1,
    outlier.shape=NA,
    font.x=font_size, font.y=font_size, font.legend=font_size, font.xtickslab=font_size, font.tickslab=font_size
    )    +
    rremove("legend") +
    rremove("ylab") +
    coord_flip() +
    ylim(c(0,40))
pdf(file.path(output.dir, "transcript_exon_number_simplified.pdf"))
print(exon_plot)
dev.off()

#combined plot
pdf(file.path(output.dir, "transcript_anno_classCode_simple_exonNum_combined.pdf"), width=14)
ggarrange(class_plot, exon_plot+rremove("y.text"), common.legend = FALSE,
          ncol = 2, nrow = 1, widths=c(1.5,1))
dev.off()

#plot mono or multi exonic in pie chart
anno_classi$exon_class <- ifelse(anno_classi$num_exons==1, "mono-exonic", "multi-exonic")
anno_classi_split <- split(anno_classi, anno_classi$class_code_simple)
exon_class <- lapply(anno_classi_split, function(x){
    as.data.frame(table(x$exon_class))
})
pal <- c("darkgray", "whitesmoke")
names(pal)<- c("multi-exonic", "mono-exonic")
labs <- paste0(exon_class[["non-chimeric (novel)"]]$Freq)
non_chimeric <- ggpie(exon_class[["non-chimeric (novel)"]], "Freq", label = labs,
        lab.pos = "in", lab.font = "black",label.size = 16, palette=pal,
        fill = "Var1", color = "white")
labs <- paste0(exon_class[["chimeric (novel)"]]$Freq)
chimeric <- ggpie(exon_class[["chimeric (novel)"]], "Freq", label = labs,
        lab.pos = "in", lab.font = "black",label.size = 16, palette=pal,
        fill = "Var1", color = "white")
labs <- paste0(exon_class[["known"]]$Freq)
known <- ggpie(exon_class[["known"]], "Freq", label = labs,
        lab.pos = "in", lab.font = "black",label.size = 16, palette=pal,
        fill = "Var1", color = "white")

pdf(file.path(output.dir, "Exon_class.pdf"), width=4)
    print(ggarrange(non_chimeric, chimeric,known,  common.legend = TRUE,
          #labels = c("non-chimeric (novel)", "chimeric (novel)", "known"),
          ncol = 1, nrow = 3
          ))
dev.off()