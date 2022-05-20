#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics_AML/comparison_gene_expression/MHCI"
dir.create(output.dir)
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
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
comp_col <- col[c(2,8,6)]
names(comp_col)<-c( "DAC_vs_DMSO", "SB939_vs_DMSO",  "DACandSB939_vs_DMSO")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
ltr_col <- c(brewer.pal(5, "Accent"), "darkgray")
names(ltr_col)<- c("LTR12","LTR12C","LTR12D","LTR12F", "no LTR12", "LTR12_")

#peptidomics list
peptides <- as.data.frame(readxl::read_excel("/omics/groups/OE0219/internal/tinat/integration/peptidomics_AML/data/AML205_218_ORF_database_I.xlsx", sheet=1))
peptides[, 23:37] <- NULL
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id

#prepare list to have a row for each accession
peptides_new <- list()
for(i in 1:length(strsplit(peptides$"Protein Accessions",";", fixed=TRUE))){
  if(length(strsplit(peptides$"Protein Accessions",";", fixed=TRUE)[[i]])==1){
    peptides_new[[i]]<- peptides[i,]
    peptides_new[[i]]$accession_new <-  peptides[i,]$"Protein Accessions"
    peptides_new[[i]]$origin<-  1
    peptides_new[[i]]$total_origins <- 1
  } else if(length(strsplit(peptides$"Protein Accessions",";", fixed=TRUE)[[i]])>1){
    peptides_new[[i]]<- peptides[rep(i,length(strsplit(peptides$"Protein Accessions",";", fixed=TRUE)[[i]])),]                    
    peptides_new[[i]]$accession_new<-  strsplit(peptides$"Protein Accessions",";", fixed=TRUE)[[i]]
    peptides_new[[i]]$origin<-  1:length(strsplit(peptides$"Protein Accessions",";", fixed=TRUE)[[i]])
    peptides_new[[i]]$total_origins <- length(strsplit(peptides$"Protein Accessions",";", fixed=TRUE)[[i]])
  }
  #if(length(strsplit(peptides_new[[i]]$accession_new,"ORF_", fixed=TRUE))>1){
  peptides_new[[i]]$transcript_id <- sapply(strsplit(peptides_new[[i]]$accession_new,"ORF_", fixed=TRUE),"[",2)
  peptides_new[[i]]$transcript_id <- sapply(strsplit( peptides_new[[i]]$transcript_id,".p", fixed=TRUE),"[",1)
  #} 
}
peptides_new <- do.call("rbind", peptides_new)
#get gene id for protein names
#library(AnnotationDbi)
#library(org.Hs.eg.db)
#mapIds(org.Hs.eg.db, keys= peptides_new$accession_new , keytype ="UNIPROT", column = "SYMBOL", multiVals = "first" )

#combine annotation with peptides
peptides_new <- dplyr::left_join(peptides_new, anno_classi, by="transcript_id")
saveRDS(peptides_new, file.path(output.dir, "peptides_list_new.rds"))
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))
write.table(peptides_new, file.path(output.dir, "peptides_list_new.tsv"), sep="\t",quote=FALSE, row.names=FALSE)

#subset peptides that originatee from our orf list
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
nrow(peptides_new_ORF)
length(unique(peptides_new_ORF$Sequence))
#select 2/3 candidates
treat_cand_seq <- peptides_new_ORF[peptides_new_ORF$"Timepoint [h]">47 ,]$Sequence
length(treat_cand_seq)
untreat_cand_seq <- peptides_new_ORF[peptides_new_ORF$"Timepoint [h]"==0 ,]$Sequence
length(untreat_cand_seq)
sel_cand <- treat_cand_seq[!treat_cand_seq %in% untreat_cand_seq]
length(sel_cand)
peptides_new_ORF_oi <- peptides_new_ORF[peptides_new_ORF$Sequence %in% sel_cand,]
nrow(peptides_new_ORF_oi)
length(unique(peptides_new_ORF_oi$Sequence))

#visualize number of transcripts peptides arise from
peptides_new_ORF_oi_sequence <- split(peptides_new_ORF_oi, peptides_new_ORF_oi$Sequence)
number_transcripts <- as.data.frame( table(sapply(peptides_new_ORF_oi_sequence, nrow))                    ) 
colnames(number_transcripts)<- c("NumberOfTranscripts", "Frequency" )
pdf(file.path(output.dir, "NumberOfTranscripts_perPeptide.pdf"))
ggbarplot(number_transcripts, x="NumberOfTranscripts", y= "Frequency",
    xlab="Number of transcripts giving rise to peptide", fill=treat_col["DAC"] )
dev.off()

#plot classification of transcripts origin
class_code_simple_stat <- as.data.frame(table(peptides_new_ORF_oi$class_code_simple))
labs <- paste0(class_code_simple_stat$Freq)
pie<- ggpie(class_code_simple_stat, "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=class_col,
                fill = "Var1", color = "black")+rremove("legend.title")
pdf(file.path(output.dir, "Pie_TranscriptClassification.pdf"), height=5, width=5)
pie
dev.off()

#plot ere anno
anno_transcript <- anno[anno$type=="transcript",]
peptides_new_ORF_oi_anno <- dplyr::left_join(peptides_new_ORF_oi, as.data.frame(anno_transcript), by="transcript_id")
peptides_new_ORF_oi_anno$ERE_TSS_anno <- NA
peptides_new_ORF_oi_anno$ERE_TSS_anno <- ifelse(peptides_new_ORF_oi_anno$dist_nearest_repeat == 0, 
    peptides_new_ORF_oi_anno$nearest_repeat_repClass, "no ERE")
peptides_new_ORF_oi_anno$ERE_TSS_anno <- ifelse(peptides_new_ORF_oi_anno$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), 
    peptides_new_ORF_oi_anno$ERE_TSS_anno , "other")
ere_stat <- as.data.frame(table(peptides_new_ORF_oi_anno$ERE_TSS_anno))
labs <- paste0(ere_stat$Freq)
pie<- ggpie(ere_stat, "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=ere_col,
                fill = "Var1", color = "black")+rremove("legend.title")
pdf(file.path(output.dir, "Pie_EREanno.pdf"), height=5, width=5)
pie
dev.off()

#plot ltr anno
peptides_new_ORF_oi_anno$LTR_TSS_anno <- NA
peptides_new_ORF_oi_anno$LTR_TSS_anno <- ifelse(peptides_new_ORF_oi_anno$dist_nearest_LTR12repeat == 0, 
    peptides_new_ORF_oi_anno$nearest_LTR12repeat_repName, "no LTR12")
ltr_stat <- as.data.frame(table(peptides_new_ORF_oi_anno$LTR_TSS_anno))
labs <- paste0(ltr_stat$Freq)
pie<- ggpie(ltr_stat, "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=ltr_col,
                fill = "Var1", color = "black")+rremove("legend.title")
pdf(file.path(output.dir, "Pie_LTRanno.pdf"), height=5, width=5)
pie
dev.off()


#Plot results in Volcano and MA plot {.tabset}
DEG_results_list_plot <- lapply(DEG_results_list, function(x){
  x <- as.data.frame(x)
  #x <- dplyr::left_join(x, anno_classi,by="transcript_id")
  x$col <- "not_significant"
  x$col <- ifelse(test = x$padj < alpha & abs(x$log2FoldChange)>lfc, yes = "significant", 
                  no ="not_significant")
  x$Label <- NA
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
    geom_point(data = DEG_results_list_plot[[i]][DEG_results_list_plot[[i]]$transcript_id %in% peptides_new_ORF_oi_anno$transcript_id,], 
        fill = "black",color = "black", alpha = 1, size=1.5) +
    ggrepel::geom_text_repel(aes(label = Label),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
    geom_vline(xintercept=c(-lfc,lfc), linetype = 2) +
    geom_hline(yintercept=-log10(alpha), linetype = 2)+
    xlab("Transcript expression change\n(log2 fold change)") + ylab("- log10(adj. P value)")
  dir.create(file.path(output.dir,i))
  ggsave(plot=Volcanos[[i]],file.path(output.dir,i,"volcano_Class_code_simple.pdf"),height = 5, width = 5, useDingbats = FALSE)
  ggsave(plot=Volcanos[[i]],file.path(output.dir,i,"volcano_Class_code_simple.png"),height = 5, width = 5, device="png")
  print(i)
}

#plot log2fold changes as boxplots
fold_change_stats <- lapply(DEG_results_list_plot[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACSB_vs_DMSO")], function(x){
    x <- x[x$transcript_id %in% peptides_new_ORF_oi_anno$transcript_id, ]
    x
})
for(i in names(fold_change_stats)){
    fold_change_stats[[i]]$comparison <- i
}
fold_change_stats_comb <- do.call("rbind",fold_change_stats)
fold_change_stats_comb$comparison <- factor(fold_change_stats_comb$comparison, levels= c(  "SB939_vs_DMSO", "DAC_vs_DMSO", "DACSB_vs_DMSO"))
pdf(file.path(output.dir,"Boxplot_log2FolChanges_peptideTranscripts.pdf"),height = 5, width = 5)
ggboxplot(fold_change_stats_comb, x="comparison", y="log2FoldChange",fill="comparison", palette= comp_col) +
        ggpubr::rotate()+
        rremove("legend") +
        rremove("ylab")
dev.off()

pdf(file.path(output.dir,"Scatter_log2FolChanges_peptideTranscripts.pdf"),height = 5, width = 5)
ggscatter(fold_change_stats_comb, x="comparison", y="log2FoldChange",fill="comparison",color="comparison", palette= comp_col) +
        geom_hline(yintercept=0, linetype = 2)+
        ggpubr::rotate()+
        rremove("legend") +
        rremove("ylab")
dev.off()

#redo boxplot but with all log2fold changes plotted besides
#plot log2fold changes as boxplots
fold_change_stats <- list()
for(i in names(DEG_results_list_plot[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACSB_vs_DMSO")])){
    fold_change_stats[[i]] <- DEG_results_list_plot[[i]][,c("transcript_id", "padj", "log2FoldChange")]
    fold_change_stats[[i]]$comparison <- i
    fold_change_stats[[i]]$candidate<-NA
    fold_change_stats[[i]]$candidate <- ifelse(fold_change_stats[[i]]$transcript_id %in%peptides_new_ORF_oi_anno$transcript_id, 
        "Yes", "No" )
}
fold_change_stats_comb <- do.call("rbind",fold_change_stats)
fold_change_stats_comb$comparison <- factor(fold_change_stats_comb$comparison, levels= c(  "SB939_vs_DMSO", "DAC_vs_DMSO", "DACSB_vs_DMSO"))
pdf(file.path(output.dir,"Boxplot_log2FolChanges_peptideTranscripts_withAllTranscripts.pdf"),height = 5, width = 5)
ggboxplot(fold_change_stats_comb, x="comparison", y="log2FoldChange",fill="candidate", legend="Peptide candidate (2/3)") +
        geom_hline(yintercept=0, linetype = 2)+
        ggpubr::rotate()+
        rremove("ylab")
dev.off()
 

# #upset plot with induced degs
# #define padj cutoff for plotting --> only up
# cutoff <- alpha
# l2fc <- lfc
# DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACSB_vs_DMSO")]
# genes2plot <- lapply(DEG_results_list_sub, function(x){
#   x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
#   x <- rownames(x)
#   x
# })
# # find non-complete elements
# ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# # remove found elements 