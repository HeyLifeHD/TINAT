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
library(rtracklayer)
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
anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
repeats<- readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))

#get transcripts of anno
anno_transcripts <- 
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

#import ORF fasta
ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")
pepteide_sub_induced <- readAAStringSet( "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")


#subset peptides that originatee from our orf list
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
#select 2/3 candidates
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]$sequences))
peptides_new_ORF_oi <- peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]
 
#match patterns
## subset orfs
peptides_oi <- AAStringSet(unique(peptides_new_ORF_oi$sequence))
names(peptides_oi)<- unique(peptides_new_ORF_oi$sequence)
#no mismatch
pattern_match_noMM <- lapply(peptides_oi, function(x){
    pattern_match <- vmatchPattern(x, pepteide_sub_induced)
    unlist(pattern_match )
})
#1 mismatch
pattern_match_1MM <- lapply(peptides_oi, function(x){
    pattern_match <- vmatchPattern(x, pepteide_sub_induced,max.mismatch=1)
    unlist(pattern_match )
})
#2 mismatch
pattern_match_2MM <- lapply(peptides_oi, function(x){
    pattern_match <- vmatchPattern(x, pepteide_sub_induced,max.mismatch=2)
    unlist(pattern_match )
})
#3 mismatch
pattern_match_3MM <- lapply(peptides_oi, function(x){
    pattern_match <- vmatchPattern(x, pepteide_sub_induced,max.mismatch=3)
    unlist(pattern_match )
})

number_of_MM <- data.frame( "Mismatches allowed: 0"=sapply(pattern_match_noMM, length),
    "Mismatches allowed: 1"=sapply(pattern_match_1MM, length),
    "Mismatches allowed: 2"=sapply(pattern_match_2MM, length),
    "Mismatches allowed: 3"=sapply(pattern_match_3MM, length)
)

number_of_MM_long <- tidyr::pivot_longer(number_of_MM, cols=colnames(number_of_MM), names_to="Number_of_mismatches_allowed", values_to="Number_of_matches_with_pepteide_sub_induced")

pdf(file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced.pdf"), height=5, width=5)
gghistogram(number_of_MM_long,fill="gray",x= "Number_of_matches_with_pepteide_sub_induced" , facet.by="Number_of_mismatches_allowed", bins=50)
dev.off()
