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
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
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
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))

#import ORF fasta
#ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")
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
    "Mismatches allowed: 2"=sapply(pattern_match_2MM, length)#,
    #"Mismatches allowed: 3"=sapply(pattern_match_3MM, length)
)

number_of_MM_long <- tidyr::pivot_longer(number_of_MM, cols=colnames(number_of_MM), names_to="Number_of_mismatches_allowed", values_to="Number_of_matches_with_pepteide_sub_induced")

pdf(file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs.pdf"), height=5, width=7)
gghistogram(number_of_MM_long[number_of_MM_long$Number_of_matches_with_pepteide_sub_induced!=0,],fill="gray",x= "Number_of_matches_with_pepteide_sub_induced" , facet.by="Number_of_mismatches_allowed", bins=50)
dev.off()

number_of_MM$peptide <- rownames(number_of_MM)
write.table(number_of_MM, file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs.txt"),row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")


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

pattern_match_df <- rbind(pattern_match_noMM_df, pattern_match_1MM_df, pattern_match_2MM_df)
#pattern_match_df <- rbind(pattern_match_0MM_df, pattern_match_1MM_df, pattern_match_2MM_df, pattern_match_3MM_df)
pattern_match_df$peptide <- sapply(strsplit(pattern_match_df$peptide, ".", fixed=TRUE), "[",1)

write.table(pattern_match_df, file.path(output.dir, "2of3DACSB_unique", "peptide_blast_pepteide_sub_induced_compared_to_inducedORFs.txt"),row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")



#scp -r heyj@172.22.62.205:/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression/2of3DACSB_unique ./

#Merge with transcript annotation table
number_of_MM <- read.delim( file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs.txt"), header=TRUE)
pattern_match_df <- read.delim(file.path(output.dir, "2of3DACSB_unique", "peptide_blast_pepteide_sub_induced_compared_to_inducedORFs.txt"),header=TRUE)
peptide_anno <- readxl::read_excel("/omics/groups/OE0219/internal/tinat/integration/peptidomics/data/220301_Peptide_Results_new.xlsx")
ORF_anno <- read.delim(file.path("/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression", "manuscript_table.tsv"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
#number_of_MM$sequences <- as.character(number_of_MM$peptide)
#number_of_MM$sequences <- gsub(" ","", number_of_MM$sequences)
#number_of_MM$peptide <- NULL
#table(number_of_MM$sequences %in% peptide_anno$sequences)
#number_of_MM <- dplyr::left_join(number_of_MM, peptide_anno)
#number_of_MM <- dplyr::left_join(number_of_MM, ORF_anno[,c("sequences", "transcript_id")])
#pattern_match_df_annonumber_of_MM <- dplyr::left_join(number_of_MM, as.data.frame(anno))
#write.table(number_of_MM, file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs_annotated.txt"),row.names = FALSE, col.names=TRUE, quote=FALSE, sep=", ")

temp <- sapply(strsplit(as.character(pattern_match_df$names), ".p", fixed=TRUE), "[", 1)
temp <- sapply(strsplit(temp, "ORF_", fixed=TRUE), "[", 2)
pattern_match_df$transcript_id <- temp

pattern_match_df_anno <- dplyr::left_join(pattern_match_df[,!colnames(pattern_match_df)%in% c("start", "end", "width")], as.data.frame(anno)[which( anno$type == "transcript"),])
write.table(pattern_match_df_anno, file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs_annotated.txt"),row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#plot ERE annotation
temp <-pattern_match_df_anno[pattern_match_df_anno$Mismatches_allowed =="0",]
dim(temp)
temp <- temp[which(temp$dist_nearest_LTR12repeat ==0),]
dim(temp)
table(temp$nearest_repeat_repName)
pattern_match_df_anno$ERE_TSS_anno <- NA
pattern_match_df_anno$ERE_TSS_anno <- ifelse(pattern_match_df_anno$dist_nearest_repeat == 0, pattern_match_df_anno$nearest_repeat_repClass, "no ERE")
pattern_match_df_anno$ERE_TSS_anno <- ifelse(pattern_match_df_anno$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), pattern_match_df_anno$ERE_TSS_anno , "other")
#split
pattern_match_df_anno_split <- split(pattern_match_df_anno, pattern_match_df_anno$Mismatches_allowed)
pattern_match_df_anno_split_stats <- lapply(pattern_match_df_anno_split, function(y){
    y <- as.data.frame(table(y$ERE_TSS_anno))
    colnames(y) <- c("repeat_class", "Freq")
    y
    })
pie <- list()
    for(j in names(pattern_match_df_anno_split_stats)){
        labs <- paste0(pattern_match_df_anno_split_stats[[j]]$Freq)
        pie[[j]] <- ggpie(pattern_match_df_anno_split_stats[[j]], "Freq", #label = labs,
                lab.pos = "in", lab.font = "black", palette=ere_col,#label.size = 16,
                fill = "repeat_class", color = "black")+rremove("legend.title")
    }
pie_eres <- ggarrange(pie[[1]], pie[[2]], pie[[3]],
    labels=paste0(names(pie), " Mismatches"),
    common.legend = TRUE,
    ncol = 1, nrow = 3
    )
pdf( file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs_annotated_Pie_transcript_class_ERE_anno.pdf"), width=5, height=10)
print(pie_eres)
dev.off()


#plot LTR annotation
pattern_match_df_anno$LTR_TSS_anno <- NA
pattern_match_df_anno$LTR_TSS_anno <- ifelse(pattern_match_df_anno$dist_nearest_LTR12repeat == 0, 
    pattern_match_df_anno$nearest_LTR12repeat_repName, "no LTR12")
#split
pattern_match_df_anno_split <- split(pattern_match_df_anno, pattern_match_df_anno$Mismatches_allowed)
pattern_match_df_anno_split_stats <- lapply(pattern_match_df_anno_split, function(y){
    y <- as.data.frame(table(y$LTR_TSS_anno))
    colnames(y) <- c("repeat_class", "Freq")
    y
    })
pie <- list()
    for(j in names(pattern_match_df_anno_split_stats)){
        labs <- paste0(pattern_match_df_anno_split_stats[[j]]$Freq)
        pie[[j]] <- ggpie(pattern_match_df_anno_split_stats[[j]], "Freq", #label = labs,
                lab.pos = "in", lab.font = "black", palette=ltr_col,#label.size = 16,
                fill = "repeat_class", color = "black")+rremove("legend.title")
    }
pie_eres <- ggarrange(pie[[1]], pie[[2]], pie[[3]],
    labels=paste0(names(pie), " Mismatches"),
    common.legend = TRUE,
    ncol = 1, nrow = 3
    )
pdf( file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs_annotated_Pie_transcript_class_LTR_anno.pdf"), width=5, height=10)
print(pie_eres)
dev.off()



#combined plotting with uniprot
pattern_match_df_anno_1 <- read.delim(file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_uniprot_annotated.txt"))
pattern_match_df_anno_2 <-  read.delim(file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_pepteide_sub_induced_compared_to_inducedORFs_annotated.txt"))

pattern_match_df_anno <- do.call(
    "rbind", list(pattern_match_df_anno_1[, c("Mismatches_allowed", "dist_nearest_repeat", "nearest_repeat_repClass","nearest_repeat_repName", "dist_nearest_LTR12repeat", "nearest_LTR12repeat_repName")], 
    pattern_match_df_anno_2[, c("Mismatches_allowed",  "dist_nearest_repeat","nearest_repeat_repClass","nearest_repeat_repName",  "dist_nearest_LTR12repeat", "nearest_LTR12repeat_repName")] ))

pattern_match_df_anno$dist_nearest_LTR12repeat <- gsub(" ", "",pattern_match_df_anno$dist_nearest_LTR12repeat)
pattern_match_df_anno$dist_nearest_repeat <- gsub(" ", "",pattern_match_df_anno$dist_nearest_repeat)

#library(dplyr)
#pattern_match_df_anno <-pattern_match_df_anno %>%
#  mutate_all(as.character)
#plot ERE annotation
temp <-pattern_match_df_anno[pattern_match_df_anno$Mismatches_allowed =="0",]
dim(temp)
temp <- temp[which(temp$dist_nearest_LTR12repeat ==0),]
dim(temp)
table(temp$nearest_repeat_repName)
pattern_match_df_anno$ERE_TSS_anno <- NA
pattern_match_df_anno$ERE_TSS_anno <- ifelse(pattern_match_df_anno$dist_nearest_repeat == 0, as.character(pattern_match_df_anno$nearest_repeat_repClass), "no ERE")
pattern_match_df_anno$ERE_TSS_anno <- ifelse(pattern_match_df_anno$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), pattern_match_df_anno$ERE_TSS_anno , "other")
#split
pattern_match_df_anno_split <- split(pattern_match_df_anno, pattern_match_df_anno$Mismatches_allowed)
pattern_match_df_anno_split_stats <- lapply(pattern_match_df_anno_split, function(y){
    y <- as.data.frame(table(y$ERE_TSS_anno))
    colnames(y) <- c("repeat_class", "Freq")
    y
    })
pie <- list()
    for(j in names(pattern_match_df_anno_split_stats)){
        labs <- paste0(pattern_match_df_anno_split_stats[[j]]$Freq)
        pie[[j]] <- ggpie(pattern_match_df_anno_split_stats[[j]], "Freq", #label = labs,
                lab.pos = "in", lab.font = "black", palette=ere_col,#label.size = 16,
                fill = "repeat_class", color = "black")+rremove("legend.title")
    }
pie_eres <- ggarrange(pie[[1]], pie[[2]], pie[[3]],
    labels=paste0(names(pie), " Mismatches"),
    common.legend = TRUE,
    ncol = 1, nrow = 3
    )
pdf( file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_combined_annotated_Pie_transcript_class_ERE_anno.pdf"), width=5, height=10)
print(pie_eres)
dev.off()


#plot LTR annotation
pattern_match_df_anno$LTR_TSS_anno <- NA
pattern_match_df_anno$LTR_TSS_anno <- ifelse(pattern_match_df_anno$dist_nearest_LTR12repeat == 0, 
    as.character(pattern_match_df_anno$nearest_LTR12repeat_repName), "no LTR12")
#split
pattern_match_df_anno_split <- split(pattern_match_df_anno, pattern_match_df_anno$Mismatches_allowed)
pattern_match_df_anno_split_stats <- lapply(pattern_match_df_anno_split, function(y){
    y <- as.data.frame(table(y$LTR_TSS_anno))
    colnames(y) <- c("repeat_class", "Freq")
    y
    })
pie <- list()
    for(j in names(pattern_match_df_anno_split_stats)){
        labs <- paste0(pattern_match_df_anno_split_stats[[j]]$Freq)
        pie[[j]] <- ggpie(pattern_match_df_anno_split_stats[[j]], "Freq", #label = labs,
                lab.pos = "in", lab.font = "black", palette=ltr_col,#label.size = 16,
                fill = "repeat_class", color = "black")+rremove("legend.title")
    }
pie_eres <- ggarrange(pie[[1]], pie[[2]], pie[[3]],
    labels=paste0(names(pie), " Mismatches"),
    common.legend = TRUE,
    ncol = 1, nrow = 3
    )
pdf( file.path(output.dir, "2of3DACSB_unique", "number_of_mismatches_combined_annotated_Pie_transcript_class_LTR_anno.pdf"), width=5, height=10)
print(pie_eres)
dev.off()





/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression/2of3DACSB_unique