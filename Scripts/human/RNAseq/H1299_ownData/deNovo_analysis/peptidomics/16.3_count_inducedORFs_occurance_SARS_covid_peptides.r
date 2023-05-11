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
library(GenomicFeaturs)

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
#anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
#anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
#anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
#anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
repeats<- readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))

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
peptides_new <- read.csv("/omics/groups/OE0219/internal/tinat/integration/peptidomics/data/sars_covid_peptides.tsv", sep="\t")

#import ORF fasta
#ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")
pepteide_sub_induced <- readAAStringSet( "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")

#subset peptides that originatee from our orf list
peptides_new_ORF <-AAStringSet(peptides_new$Sequence )
names(peptides_new_ORF) <- peptides_new$Peptide.ID

#match patterns
## subset orfs
peptides_oi <- peptides_new_ORF

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

dir.create(file.path(output.dir, "SARS_Covid_peptide_match"))
pdf(file.path(output.dir, "SARS_Covid_peptide_match", "number_of_mismatches_SARS_Covid_peptide_compared_to_inducedORFs.pdf"), height=5, width=7)
gghistogram(number_of_MM_long[number_of_MM_long$Number_of_matches_with_pepteide_sub_induced!=0,],fill="gray",x= "Number_of_matches_with_pepteide_sub_induced" , facet.by="Number_of_mismatches_allowed", bins=50)
dev.off()

number_of_MM$peptide <- rownames(number_of_MM)
write.table(number_of_MM, file.path(output.dir, "SARS_Covid_peptide_match", "number_of_mismatches_SARS_Covid_peptide_compared_to_inducedORFs.txt"),row.names = FALSE, col.names=TRUE, quote=FALSE, sep=", ")


#generate detailed table
pattern_match_noMM_df <- do.call("rbind",lapply(pattern_match_noMM, function(x){
    x <- as.data.frame(x)
    x
}))

pattern_match_noMM_df$peptide <- rownames(pattern_match_noMM_df)
pattern_match_noMM_df$Mismatches_allowed <- "0"

pattern_match_1MM_df <- do.call("rbind",lapply(pattern_match_1MM, function(x){
    x <- as.data.frame(x)
    x
}))
pattern_match_1MM_df$peptide <- rownames(pattern_match_1MM_df)
pattern_match_1MM_df$Mismatches_allowed <- "1"

pattern_match_2MM_df <- do.call("rbind",lapply(pattern_match_2MM, function(x){
    x <- as.data.frame(x)
    x
}))
pattern_match_2MM_df$peptide <- rownames(pattern_match_2MM_df)
pattern_match_2MM_df$Mismatches_allowed <- "2"

pattern_match_3MM_df <- do.call("rbind",lapply(pattern_match_3MM, function(x){
    x <- as.data.frame(x)
    x
}))
pattern_match_3MM_df$peptide <- rownames(pattern_match_3MM_df)
pattern_match_3MM_df$Mismatches_allowed <- "3"

#pattern_match_df <- rbind(pattern_match_noMM_df, pattern_match_1MM_df, pattern_match_2MM_df)
pattern_match_df <- rbind(pattern_match_noMM_df, pattern_match_1MM_df, pattern_match_2MM_df, pattern_match_3MM_df)
pattern_match_df$peptide <- sapply(strsplit(pattern_match_df$peptide, ".", fixed=TRUE), "[",1)

write.table(pattern_match_df, file.path(output.dir, "SARS_Covid_peptide_match", "peptide_blast_pepteide_SARS_Covid_peptide_compared_to_inducedORFs.txt"),row.names = FALSE, col.names=TRUE, quote=FALSE, sep=", ")



/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression/SARS_Covid_peptide_match ./

scp -r heyj@172.22.62.205:/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression/SARS_Covid_peptide_match ./

