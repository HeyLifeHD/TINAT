#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_analysis"
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
counts <- read.table("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/transcripts/transcript_count_matrix.csv", sep=",",  header=TRUE)
rownames(counts)<- counts$transcript_id
counts$transcript_id <- NULL
colnames(counts)<-  sapply(strsplit( colnames(counts),".sorted",fixed=TRUE), "[",1)
colnames(counts) <- gsub(".", "_", colnames(counts),fixed=TRUE)

#subset counts
#subset failed samples
counts <- counts[,! colnames(counts) %in% c("AS_811024_LR_62189", "AS_811004_LR_62189")]
#create annotation sheet
files <- list.files("/omics/odcf/project/OE0219/tinat/sequencing/rna_sequencing/view-by-pid/OE0219_TINAT_B16-F10/", pattern="fastq.gz",full.names= TRUE, recursive=TRUE)
files <- files[grep("_R1", files)]
files <- sapply(strsplit( files,"//",fixed=TRUE), "[",2)
file_names <- sapply(strsplit( files,"/",fixed=TRUE), "[",5)
file_names <- sapply(strsplit( file_names,"_R1",fixed=TRUE), "[",1)
files <- sapply(strsplit( files,"/",fixed=TRUE), "[",1)
names(files)<- file_names
names(files)
names(files) <- gsub("-", "_", names(files),fixed=TRUE)
files_sub <- files[colnames(counts)]

treatment= sapply(strsplit( files_sub,"-",fixed=TRUE), "[",1)
treatment <- toupper(treatment)
length(treatment)
replicate= sapply(strsplit( files_sub,"-",fixed=TRUE), "[",2)
replicate <- paste0("rep_", replicate)
length(replicate)
Sample_ID=as.character(paste0(treatment,paste0("_Rep", replicate)))
sample_anno <- data.frame(Sample_ID=as.character(Sample_ID), treatment=treatment, replicate=as.factor(replicate))
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