#plot feature counts of repeats -SRP230446 - paired end
#libraries
library(rtracklayer)
library(data.table)
library(DESeq2)

#folders
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220622_AMLPatient_FeatureCountRepeat_analysis_revised/"
dir.create(base.dir)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)

#load counts
counts <- fread("/omics/groups/OE0219/internal/tinat/Cellline_panel/220608_AMLPatient_FeatureCountsRepeats/repeats_featureCounts.txt")
anno <- counts[,1:6]
counts <- counts[,7:ncol(counts)]
#rename counts sample names
#create annotation sheet

#create annotation sheet
files <- colnames(counts)
files <- sapply(strsplit(files, "/"),"[",11)
files <- sapply(strsplit(files, ".", fixed=TRUE),"[",1)

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
colnames(counts)==rownames(sample_anno)
#subset qc failed samples
counts <- counts[,!colnames(counts) %in% c("A14001_c1d8", "A26004_c1d8", "A03005_c1d8"),with=FALSE ]
dim(counts)
sample_anno <- sample_anno[rownames(sample_anno)%in%colnames(counts), ]
dim(sample_anno)

#counts <- counts[,rownames(sample_anno)]
#import_repeats
repeats<-import.gff("/omics/groups/OE0219/internal/repguide/data/repeats/repeats_hg19.gtf.gz")#readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))
#import anno
counts_repeats <- cbind(as.data.frame(counts), repeats[, c("gene_id", "transcript_id")])
#count repeat families and classes
counts_sum_repFamily <- aggregate(counts_repeats[,1:(ncol(counts_repeats)-7)], by=list(repFamily=counts_repeats$gene_id), FUN=sum)
counts_sum_repName <- aggregate(cbind(counts_repeats[,1:(ncol(counts_repeats)-7)]), by=list(repFamily=counts_repeats$transcript_id), FUN=sum)
rownames(counts_sum_repFamily) <- counts_sum_repFamily$repFamily
counts_sum_repFamily[,1]<- NULL
counts_sum_repFamily <- counts_sum_repFamily[,rownames(sample_anno)]

#greate dds
dds_repeats <- DESeqDataSetFromMatrix(countData = counts_sum_repFamily, 
                              colData = sample_anno, 
                              design =  ~  Patient_ID + group )

#Estimating size factors
dds_repeats <- estimateSizeFactors(dds_repeats)

#Running the differential expression 
dds_repeats <- DESeq(dds_repeats)
saveRDS(dds_repeats, file.path(results.dir, "dds_repFamily_group.rds"))
