#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220622_AMLPatient_deNovoCellline_quantification_analysis_revised/"
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
counts <- read.table("/omics/groups/OE0219/internal/tinat/Cellline_panel/220607_AMLPatient_deNovoCellline_quantification/transcripts/transcript_count_matrix.csv", sep=",",  header=TRUE)
rownames(counts)<- counts$transcript_id
counts$transcript_id <- NULL
dim(counts)

#create annotation sheet
files <- colnames(counts)
files_sub <-  strsplit(files, "_", fixed=TRUE)
#SampleID <- sapply(files_sub, function(x)paste0(x[2],"_", x[3]))
treatment <-  sapply(files_sub, function(x)paste0(x[3]))
treatment <- gsub(".sorted", "",treatment)
patientID <-  sapply(files_sub, function(x)paste0(x[2]))
SampleID <- paste0(patientID,"_", treatment)

sample_anno <- data.frame(Sample_ID=as.character(SampleID), Patient_ID=patientID,treatment=treatment)
rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)
sample_anno$group <- ifelse(sample_anno$treatment %in% c("scr", "c1d1"), "pre_treatment", "post_treatment")
#rename counts
colnames(counts) <- rownames(sample_anno)
counts <- counts[,rownames(sample_anno)]
#subset qc failed samples
counts <- counts[,!colnames(counts) %in% c("A14001_c1d8", "A26004_c1d8", "A03005_c1d8") ]
dim(counts)
sample_anno <- sample_anno[rownames(sample_anno)%in%colnames(counts), ]
dim(sample_anno)
#create dds
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_anno, 
                              design = ~ Patient_ID + group )

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


# sample_anno$library <- ifelse(sample_anno$cellline=="H1299_H2", "run1", "run2")
# dds <- DESeqDataSetFromMatrix(countData = counts, 
#                               colData = sample_anno, 
#                               design = ~library+ cellline + treatment )


