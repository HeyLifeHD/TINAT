#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)
library(Biostrings)

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"
dir.create(output.dir)
#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))

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
ltr_col <- brewer.pal(5, "Accent")
names(ltr_col)<- c("LTR12","LTR12C","LTR12D","LTR12F", "no LTR12")
#peptidomics list
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))

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

#import ORF fasta
ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs.pep")


#get name of orfs and width as dataframe
ORF_stats <- data.frame(ORF_names=as.character(names(ORFs)), width=width(ORFs))
#temp <- sapply(strsplit(as.character(ORF_stats$ORF_names), " type",fixed=TRUE),"[",1) 
#temp <- paste0("ORF_",sapply(strsplit(temp, " type",fixed=TRUE), "[",1))
ORF_stats$ORF_name_short <- sapply(strsplit(as.character(ORF_stats$ORF_names), " type",fixed=TRUE),"[",1) 
ORF_stats$transcript_id <- sapply(strsplit(as.character(ORF_stats$ORF_name_short), ".p",fixed=TRUE),"[",1) 

#combine
ORF_stats <- dplyr::left_join(ORF_stats, anno_classi, by="transcript_id")

#spliut based on class
ORF_stats_split <- split(ORF_stats, ORF_stats$class_code_simple)

for(i in names(ORF_stats_split)){
#order based on width
 ORF_stats_split[[i]] <-  ORF_stats_split[[i]][order( ORF_stats_split[[i]]$width, decreasing=TRUE),]
 ORF_stats_split[[i]]$order <- 1:nrow( ORF_stats_split[[i]])

#plot as rainfall plot
plot <- ggplot( ORF_stats_split[[i]], aes(x=width, y=order))+
    theme_pubr()+
    geom_point(data =  ORF_stats_split[[i]], 
        fill = "gray",color = "gray", alpha = 0.5, size=.1) +
    geom_vline(xintercept=mean( ORF_stats_split[[i]]$width), linetype = 2, color="gray") +
    annotate(geom="text", x=mean( ORF_stats_split[[i]]$width)+10, y=median( ORF_stats_split[[i]]$order), 
        label=paste0("Mean ORF length: ", round(mean( ORF_stats_split[[i]]$width)), " AAs"), color="gray") +
    scale_x_log10()+
    xlab("Length of ORF [AA]") + 
    ylab("Order by length")

pdf(file.path(output.dir, paste0("ORF_length_allORFs_mean_",i, ".pdf")), height=5, width=5)
print(ggrastr::rasterize(plot, layers="Point", dpi=100))
dev.off()
print(i)
}