#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
dir.create(base.dir)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)

#load data
counts <- read.table("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification/transcripts/transcript_count_matrix.csv", sep=",",  header=TRUE)
rownames(counts)<- counts$transcript_id
counts$transcript_id <- NULL
#create annotation sheet
treatment= c(rep(c("DMSO", "DAC", "SB939", "DACandSB939"),3))
length(treatment)
replicate=c(rep(1, 4), rep(2, 4), rep(3, 4) )
length(replicate)
Sample_ID=as.character(paste0(treatment,paste0("_Rep", replicate)))
sample_anno <- data.frame(Sample_ID=as.character(Sample_ID), treatment=treatment, replicate=as.factor(replicate))

temp <- sample_anno
temp$Sample_ID <- sapply(strsplit(colnames(counts), ".sorted"), "[", 1)
temp$Sample_ID <- gsub(".", "-", temp$Sample_ID, fixed=TRUE)
temp <- temp[order(temp$Sample_ID),]
write.table(temp, file.path(results.dir, "sample_anno.csv"), sep=";", row.names=FALSE, col.names=TRUE, quote=FALSE)

rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)

#rename counts
colnames(counts) <-sample_anno$Sample_ID

#create dds
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_anno, 
                              design = ~ treatment )

#Estimating size factors
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=FALSE) >= 1 ) >= 1
table(idx)
#dds <- dds[idx,]
#dim(dds)

#Running the DESEQ
dds <- DESeq(dds)
saveRDS(dds, file.path(results.dir,"dds.rds"))