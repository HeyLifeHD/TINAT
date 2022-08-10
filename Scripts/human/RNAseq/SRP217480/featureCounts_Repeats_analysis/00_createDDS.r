#plot feature counts of repeats -SRP230446 - paired end
#libraries
library(rtracklayer)
library(data.table)
library(DESeq2)

#folders
base.dir<- "/omics/groups/OE0219/internal/tinat/211007_SRP217480_featureCount_Repeats/"
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
counts <- fread("/omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/repeats_featureCounts.txt")
anno <- counts[,1:6]
counts <- counts[,7:ncol(counts)]

colnames(counts)<- sapply(strsplit(colnames(counts), "/", fixed=TRUE),"[", 10)
colnames(counts)<- sapply(strsplit(colnames(counts), "_", fixed=TRUE),"[", 1)

#rename counts sample names
#create annotation sheet
sample_anno <- as.data.frame(fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP217480/sample_sheet_new.txt"))
sample_anno <- sample_anno[!duplicated(sample_anno),]
cellline <- sapply(strsplit(sample_anno$organism_taxid, ".", fixed=TRUE),"[", 1)
concentration <- sapply(strsplit(sample_anno$organism_taxid, ".", fixed=TRUE),"[", 2)
treatment <- sapply(strsplit(sample_anno$organism_taxid, ".", fixed=TRUE),"[", 3)
day <- sapply(strsplit(sample_anno$organism_taxid, ".", fixed=TRUE),"[", 4)
replicate <- sapply(strsplit(sample_anno$organism_taxid, ".", fixed=TRUE),"[", 5)

sample_anno <- data.frame(Sample_ID=as.character(sample_anno$sample_id), cellline=cellline, 
    treatment=treatment,concentration=concentration,day=day, replicate=as.factor(replicate), 
    grouping=paste0(cellline, "_", treatment,"_", concentration,"_",day))
rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)
#rename counts
sample_anno <- sample_anno[colnames(counts), ]


#import_repeats
repeats<-import.gff("/omics/groups/OE0219/internal/repguide/data/repeats/repeats_hg19.gtf.gz")#readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))
#import anno
counts_repeats <- cbind(as.data.frame(counts), repeats[, c("gene_id", "transcript_id")])
#count repeat families and classes
counts_sum_repFamily <- aggregate(counts_repeats[,1:(ncol(counts_repeats)-7)], by=list(repFamily=counts_repeats$gene_id), FUN=sum)
#counts_sum_repName <- aggregate(cbind(counts_repeats[,1:(ncol(counts_repeats)-7)]), by=list(repFamily=counts_repeats$transcript_id), FUN=sum)
rownames(counts_sum_repFamily) <- counts_sum_repFamily$repFamily
counts_sum_repFamily[,1]<- NULL

#greate dds
dds_repeats <- DESeqDataSetFromMatrix(countData = counts_sum_repFamily, 
                              colData = sample_anno, 
                              design = ~  grouping)

#Estimating size factors
dds_repeats <- estimateSizeFactors(dds_repeats)

#Running the differential expression 
dds_repeats <- DESeq(dds_repeats)
saveRDS(dds_repeats, file.path(results.dir, "dds_repFamily_group.rds"))
