#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/Encode/220128_deNovo_quantification_H1299ref_analysis/"
dir.create(base.dir,recursive=TRUE)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)


#load data
#for reverse
counts1 <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_reverse/transcripts/transcript_count_matrix.csv"))
#for unstranded
counts2 <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded/transcripts/transcript_count_matrix.csv"))

#combine count tables
# if(rownames(counts1)==rownames(counts2)){
#     rownames(counts2)<- counts2$transcript_id
#     rownames(counts1)<- counts1$transcript_id
#     counts2$transcript_id <- NULL
#     counts1$transcript_id <- NULL
#     counts <- cbind(counts1, counts2)
# } else {
    counts <- dplyr::inner_join(counts1, counts2, by="transcript_id")
    rownames(counts)<- counts$transcript_id
    counts$transcript_id <- NULL
#}

#redefine couts names
colnames(counts) <- gsub(".sorted_transcripts","", colnames(counts)) 

#replace 0 with na
lapply(counts, function(x){
    table(is.na(x))
})
counts[is.na(counts)]<- 0

#check counts
all.zero <- apply(counts, 2, function(x) count(x==0))
#remove counts with all zero
counts <- counts[,!colnames(counts) %in% names(all.zero[all.zero ==nrow(counts)])]

#create annotation sheet
meta <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/Encode/data_new","meta_new.rds"))
meta$sra <- meta$dbxrefs
meta$sra <- gsub("SRA:", "", meta$sra)
meta<- meta[meta$"Paired end"==1,]
rownames(meta)<- meta$sra

#check meta for quality and subset counts
#meta_sub <- meta[meta$"Audit NOT_COMPLIANT" =="",]
#nrow(meta_sub)/nrow(meta)
#meta_sub <- meta[meta$"Audit NOT_COMPLIANT" =="",]
# nrow(meta_sub)/nrow(meta)
# counts_sub <- counts[,colnames(counts)%in% rownames(meta_sub)]
# meta_sub <- meta_sub[colnames(counts_sub),]
counts_sub <- counts[,colnames(counts)%in% rownames(meta)]
meta_sub <- meta[colnames(counts_sub),]

#prepare sample_anno
sample_anno <- meta_sub[,c("sra", "Assay", "Biosample term name", "Library made from","Library strand specific",
     "Biological replicate(s)", "Technical replicate(s)", "Lab")]
colnames(sample_anno)<- c("sra", "assay", "tissue", "input", "strandedness","biological_replicate", 
    "technical_replicate", "lab")

#####size estimation failes --> add pseudocounts
#counts_sub <- counts+1


#create dds
dds <- DESeqDataSetFromMatrix(countData = counts_sub, 
                              colData = sample_anno, 
                              design = ~tissue )


#Estimating size factors
#for Genotype comparison
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
table(idx)
#dds <- dds[idx,]
#dim(dds)

#Running the DESEQ
#dds <- DESeq(dds, parallel=T^RUE, BPPARAM=MulticoreParam(4))
saveRDS(dds, file.path(results.dir,"dds.rds"))

