#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/210913_SRP073103_processing_deNovo_analysis/T98G"
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
counts <- read.table("/omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/T98G/transcripts/transcript_count_matrix.csv", sep=",",  header=TRUE)
rownames(counts)<- counts$transcript_id
counts$transcript_id <- NULL
#create annotation sheet
treatment= c(rep("DMSO",2), rep("DAC",2))
length(treatment)
replicate=c(rep(c(1,2), 2))
length(replicate)
Sample_ID=as.character(paste0(treatment,paste0("_Rep", replicate)))
sample_anno <- data.frame(Sample_ID=as.character(Sample_ID), treatment=treatment, replicate=as.factor(replicate))
rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)
#rename counts
colnames(counts) <-sample_anno$Sample_ID

#create dds
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_anno, 
                              design = ~replicate+ treatment )

#Estimating size factors
#for Genotype comparison
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
table(idx)
#dds <- dds[idx,]
#dim(dds)

#Running the DESEQ
dds <- DESeq(dds)
saveRDS(dds, file.path(results.dir,"dds.rds"))