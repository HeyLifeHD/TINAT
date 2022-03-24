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
temp <-  readRDS("/omics/groups/OE0219/internal/tinat/integration/peptidomics/data/TINAT_occupancy_matrices_revision.RDS" )
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

#combine with gene expression
DEG_results_list_sub <- lapply(DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")], function(x){
    x <- x[rownames(DEG_results_list[[1]]),]
    x[,c("padj", "log2FoldChange","nearest_repeat_repName","dist_nearest_repeat")]
})
gene_expr <- do.call("cbind",DEG_results_list_sub)
gene_expr$transcript_id <- rownames(gene_expr)
peptides_new <- dplyr::left_join(peptides_new, gene_expr, by="transcript_id")
saveRDS(peptides_new, file.path(output.dir, "peptides_list_new.rds"))
write.table(peptides_new, file.path(output.dir, "peptides_list_new.tsv"), sep="\t",quote=FALSE, row.names=FALSE)


#subset peptides that originatee from our orf list
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
#select 2/3 candidates
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB>67 & peptides_new_ORF$DMSO ==0,]$sequences))
peptides_new_ORF_oi <- peptides_new_ORF[peptides_new_ORF$DAC_SB>67 & peptides_new_ORF$DMSO ==0,]


peptides_new_ORF_oi_split <- split(peptides_new_ORF_oi, peptides_new_ORF_oi$sequences)
peptides_new_ORF_oi_split_sub <- lapply(peptides_new_ORF_oi_split, function(x){
    x <- x[order(x$DACandSB939_vs_DMSO.log2FoldChange, decreasing=TRUE),]
    x[1,]
})
peptides_new_ORF_oi_split_sub <- do.call("rbind",peptides_new_ORF_oi_split_sub)
write.table(peptides_new_ORF_oi_split_sub, file.path(output.dir, "manuscript_table.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
