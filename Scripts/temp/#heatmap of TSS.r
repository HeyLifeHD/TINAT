#overlap of tsss
#libraries
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(peakSeason)

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"
input.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/6_tracks_normalized"
#load data
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#Set specifications 
cutoff <- 0.01 #set FDR cutoff
l2fc <- 2##set logfold2 cutoff

#annotations
#anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
#anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
repeats <- fread("/omics/groups/OE0219/internal/repguide/data/repeats/hg19_repeats.txt")
responsive_TSS <- readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/responsive_TSSs.RDS")

#subset tss from david
responsive_TSS <- responsive_TSS[-2]
names(responsive_TSS) <- names(DEG_results_list_sub_TSS_novel)

#merge deg with anno class
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id 
#combined plot
anno_transcript <- anno[anno$type == "transcript",]
DEG_results_list_sub_TSS_novel <- lapply(DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")], function(x){
        x <- dplyr::left_join(x, anno_classi, by="transcript_id")
        x <- x[x$class_code_simple!="known",]
        x <- resize(anno_transcript[anno_transcript$transcript_id %in% x[which(x$padj < alpha & x$log2FoldChange>lfc),]$transcript_id,],1,ignore.strand=FALSE)
        #x$peak <- NA
        #x$peak <- paste0(1:nrow(x), "region")
        x
}) 

#look at overlap
subset_TSS <- list()
for(i in names(DEG_results_list_sub_TSS_novel)){
    dir.create(file.path(PostDE.dir,i, "TINAT_compare" ))
    temp <- responsive_TSS[[i]]
    temp <- unique(temp)
    temp <- temp+100
    
    temp2 <- DEG_results_list_sub_TSS_novel[[i]]
    mcols(temp2)<- NULL
    temp2 <- unique(temp2)
    temp2 <- temp2+100

    print(length(subset_TSS[[i]] ))
    pdf(file.path(PostDE.dir,i, "TINAT_compare", "TSS_comparison_expand100bp.pdf"))
    print(ChIPpeakAnno::makeVennDiagram(list(RNAseq=temp2, CAGE=temp), 
        NameOfPeaks=c("RNAseq", "CAGE"), by="region", ignore.strand=TRUE))
    dev.off()
}