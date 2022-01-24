#sense antisense pairs
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(knitr)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)
library(clusterProfiler)
library(readxl)
library(rtracklayer)


#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
dds <-readRDS(file.path(results.dir, "dds.rds"))
vst <- readRDS(file.path(results.dir, "vst.rds"))
anno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
#Take a look at design and annotation of samples
design(dds)
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

#subset transcripts
anno_transcripts <- anno[anno$type =="transcript",]

#focus on + strand
anno_transcripts_minus <- anno_transcripts[strand(anno_transcripts)=="-",]
#find overlaps
ol <- findOverlaps(anno_transcripts[strand(anno_transcripts)=="+",], anno_transcripts_minus, ignore.strand=TRUE)
ol_df <- as.data.frame(ol)
#split overlaps
ol_df_split <- split(ol_df, ol_df$queryHits)
names(ol_df_split) <- anno_transcripts[strand(anno_transcripts)=="+",][unique(ol_df$queryHits),]$transcript_id
#merge with info
anno_transcripts_pairs_plus <- list()
for(i in names(ol_df_split)){
    anno_transcripts_pairs_plus[[i]] <- anno_transcripts_minus[ol_df_split[[i]]$subjectHits,]
    print(which(names(ol_df_split)==i)/length( names(ol_df_split)))
}


#focus on - strand
anno_transcripts_plus <- anno_transcripts[strand(anno_transcripts)=="+",]
#find overlaps
ol <- findOverlaps(anno_transcripts[strand(anno_transcripts)=="-",], anno_transcripts_plus, ignore.strand=TRUE)
ol_df <- as.data.frame(ol)
#split overlaps
ol_df_split <- split(ol_df, ol_df$queryHits)
names(ol_df_split) <- anno_transcripts[strand(anno_transcripts)=="-",][unique(ol_df$queryHits),]$transcript_id
#merge with info
anno_transcripts_pairs_minus <- list()
for(i in names(ol_df_split)){
    anno_transcripts_pairs_minus[[i]] <- anno_transcripts_plus[ol_df_split[[i]]$subjectHits,]
    print(which(names(ol_df_split)==i)/length( names(ol_df_split)))
}
#merge
anno_transcripts_pairs<- c(anno_transcripts_pairs_plus, anno_transcripts_pairs_minus)

#merge with deg analysis
anno_transcripts_pairs_merged <- parallel::mclapply(anno_transcripts_pairs, function(x){
    x <- dplyr::left_join(as.data.frame(x),DEG_results_list$DACandSB939_vs_DMSO, by="transcript_id")
    x
}, mc.cores=10)

#check for expression patterns
anno_transcripts_pairs_merged_bothInfo <- list()
for(i in names(anno_transcripts_pairs_merged)){
    anno_transcripts_pairs_merged_bothInfo[[i]]<- list()
    anno_transcripts_pairs_merged_bothInfo[[i]]$query <- DEG_results_list$DACandSB939_vs_DMSO[DEG_results_list$DACandSB939_vs_DMSO$transcript_id ==i, ]
    anno_transcripts_pairs_merged_bothInfo[[i]]$subject <- anno_transcripts_pairs_merged[[i]]
    print(which(names(anno_transcripts_pairs_merged)==i)/length( names(anno_transcripts_pairs_merged)))
}
#subset query based on lfc and alpha 
index_querySub <- sapply(anno_transcripts_pairs_merged_bothInfo, function(x){
    x <- ifelse(x[["query"]]$padj < alpha & abs(x[["query"]]$log2FoldChange)>lfc, TRUE, FALSE)
    x
})
anno_transcripts_pairs_merged_bothInfo_querySub <- anno_transcripts_pairs_merged_bothInfo[index_querySub]
length(anno_transcripts_pairs_merged_bothInfo_querySub)/length(anno_transcripts_pairs_merged_bothInfo)

#subset subject based on opposite trend 
index_subjectSub <- lapply(anno_transcripts_pairs_merged_bothInfo_querySub, function(x){
    x <- if(x[["query"]]$log2FoldChange > 0){
        x <- ifelse(x[["subject"]]$padj < alpha & x[["subject"]]$log2FoldChange < (-lfc), TRUE, FALSE)
        }else if(x[["query"]]$log2FoldChange < 0){
        x <- ifelse(x[["subject"]]$padj < alpha & x[["subject"]]$log2FoldChange>lfc, TRUE, FALSE)
    }
    x
})
anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub <- list()
for(i in names(index_subjectSub)){
    if(sum(index_subjectSub[[i]]>=1)){
            anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub[[i]] <- list()
            anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub[[i]]$query<- anno_transcripts_pairs_merged_bothInfo_querySub[[i]][["query"]]
            anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub[[i]]$subject <- anno_transcripts_pairs_merged_bothInfo_querySub[[i]][["subject"]][index_subjectSub[[i]],]

    }
}
length(anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub)
saveRDS(anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub, file.path(PostDE.dir, "anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub.rds"))
#create simple list
df <- list()
for(i in names(anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub)){
    df[[i]] <- data.frame(query=i, subject=anno_transcripts_pairs_merged_bothInfo_querySub_subjectSub[[i]][["subject"]]$transcript_id)
    
}
df <- do.call("rbind",df)
write.table(df, file.path(PostDE.dir, "querySubject_pairs.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")