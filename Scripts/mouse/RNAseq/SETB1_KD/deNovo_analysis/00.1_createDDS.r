#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220811_RNAseqSETB1KD_deNovoB16DACSB_analysis"
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
counts <- read.table("/omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/transcripts/transcript_count_matrix.csv", sep=",",  header=TRUE)
rownames(counts)<- counts$transcript_id
counts$transcript_id <- NULL

#create anno
treatment= sapply(strsplit( colnames(counts),"_",fixed=TRUE), "[",2)
treatment <- gsub("SETB1", "SETB1_KD",treatment)
replicate <- sapply(strsplit( colnames(counts),"rep",fixed=TRUE), "[",2)
replicate <- sapply(strsplit( replicate,"_1",fixed=TRUE), "[",1)
replicate <- paste0("rep_", replicate)
gRNA <-  sapply(strsplit( colnames(counts),"sg",fixed=TRUE), "[",2)
gRNA <-  sapply(strsplit( gRNA,"_rep",fixed=TRUE), "[",1)
gRNA <- paste0("sg_", gRNA)
sample_id <- paste0(treatment, "_", gRNA, "_", replicate)

sample_anno <- data.frame(Sample_ID=as.character(sample_id), treatment=treatment,gRNA=gRNA, replicate=as.factor(replicate))
write.table(sample_anno, file.path(results.dir, "sample_anno.csv"), sep=";", row.names=FALSE, col.names=TRUE, quote=FALSE)
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