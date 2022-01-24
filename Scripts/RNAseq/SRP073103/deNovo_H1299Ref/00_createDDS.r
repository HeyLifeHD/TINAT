#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/211007_SRP073103_deNovo_H1299_analysis"
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
counts <- read.table("/omics/groups/OE0219/internal/tinat/210929_SRP073103_deNovo_quantification_H1299/transcripts/transcript_count_matrix.csv", sep=",",  header=TRUE)
rownames(counts)<- counts$transcript_id
counts$transcript_id <- NULL
#create annotation sheet
treatment= c(rep(c(rep("DMSO",2), rep("DAC",2)),3))
length(treatment)
replicate=c(rep(rep(c(1,2), 2),3))
length(replicate)
cellline = c(rep("LNT",4), rep("T98G",4), rep("U87",4))
length(cellline)
Sample_ID=as.character(paste0(cellline, "_",treatment,paste0("_Rep", replicate)))
grouping= paste0(paste0(cellline, "_",treatment))
sample_anno <- data.frame(Sample_ID=as.character(Sample_ID),grouping=grouping, treatment=treatment, replicate=as.factor(replicate))
rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)
#rename counts
colnames(counts) <-sample_anno$Sample_ID

#create dds
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_anno, 
                              design = ~grouping )

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