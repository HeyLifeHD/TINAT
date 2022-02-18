#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)

#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#function
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}
#DESeq2 Analysis
dds <-readRDS(file.path(results.dir, "dds.rds"))

#get annotation
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")

#get counts
counts_raw <-  counts(dds, normalized=FALSE, replaced=FALSE) 

#get transcript length
anno_exons <- anno[anno$type=="exon",]
anno_exons_split <- split( anno_exons,anno_exons$transcript_id)
exons_length <- lapply(anno_exons_split, function(x){
    x <- width(x)
    x <- sum(x)
})
exons_length <- unlist(exons_length)#
saveRDS(exons_length, file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4", "exon_length.rds"))
exons_length <-readRDS(file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4", "exon_length.rds"))
#calculate tpms
tpm <- counts_to_tpm(counts=counts_raw, featureLength=exons_length[rownames(counts_raw)], meanFragmentLength=rep(1, ncol(counts_raw)))
saveRDS(tpm, file.path(results.dir, "tpm_from_counts.rds"))

#summarize over mean tissue
sample_anno <- colData(dds)
sample_anno_split <- split(sample_anno, sample_anno$treatment)
tpm_mean <- sapply(sample_anno_split, function(x){
    if(length(x$sra)==1){
        x <- tpm[, rownames(x)]
    } else{
        x <- rowMeans(tpm[,rownames(x)])
    }
    x
})
saveRDS(tpm_mean, file.path(results.dir, "tpm_meanTreatment_from_counts.rds"))
