#libraries
library(DESeq2)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification_analysis/"
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
counts <- read.table("/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification/transcripts/transcript_count_matrix.csv", sep=",",  header=TRUE)
rownames(counts)<- counts$transcript_id
counts$transcript_id <- NULL
dim(counts)
#create annotation sheet
files <- list.files("/omics/odcf/project/OE0219/tinat/sequencing/rna_sequencing/view-by-pid/", 
    pattern="1.fastq.gz$", full.names=TRUE, recursive=TRUE)
files_sub <-  files[sapply(strsplit(files, "/", fixed=TRUE), function(x) {
    x[14] %in% c("run210622_VH00211_82_AAACN3VHV", "run220202_VH00211_134_AAACMNTHV",  "run220217_VH00693_26_AAACLW7HV")
    }
    )]
files_sep <- strsplit(files_sub, "/", fixed=TRUE)
treatment <- toupper(sapply(strsplit(sapply(files_sep, "[",12),"-", fixed=TRUE ),"[",1))
replicate <- paste0("rep_",sapply(strsplit(sapply(files_sep, "[",12),"-", fixed=TRUE ),"[",2))
cellline <- sapply(strsplit(sapply(files_sep, "[",11),"TINAT_", fixed=TRUE ),"[",2)
SampleID <- sapply(strsplit(sapply(files_sep, "[",16),"_R", fixed=TRUE ),"[",1)
SampleID <- gsub("-",".", SampleID)
sample_anno <- data.frame(Sample_ID=as.character(SampleID), cellline=cellline,treatment=treatment, replicate=as.factor(replicate))
rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)
#rename counts
colnames(counts) <-sapply(strsplit(colnames(counts),".sorted_"),"[",1)
counts <- counts[,rownames(sample_anno)]

#create dds
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_anno, 
                              design = ~ cellline + treatment )

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