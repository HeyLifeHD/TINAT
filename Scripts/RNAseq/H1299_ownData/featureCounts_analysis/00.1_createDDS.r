#plot feature counts of repeats -SRP230446 - paired end
#libraries
library(rtracklayer)
library(data.table)
library(DESeq2)

#folders
base.dir<- "/omics/groups/OE0219/internal/tinat/210802_shortRead_FeatureCountsRepeats_analysis/"
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
counts <- fread("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/repeats_featureCounts.txt")
anno <- counts[,1:6]
counts <- counts[,7:ncol(counts)]
#rename counts sample names
#create annotation sheet
treatment= c(rep(c("DMSO", "DAC", "SB939", "DACandSB939"),3))
length(treatment)
replicate=c(rep(1, 4), rep(2, 4), rep(3, 4) )
length(replicate)
Sample_ID=as.character(paste0(treatment,paste0("_Rep", replicate)))
sample_anno <- data.frame(Sample_ID=as.character(Sample_ID), treatment=treatment, replicate=as.factor(replicate))
rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)
#rename counts
colnames(counts) <-as.character(sample_anno$Sample_ID)


#import_repeats
repeats<-import.gff("/icgc/dkfzlsdf/analysis/C010/repguide/data/repeats/repeats_hg19.gtf.gz")#readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))
#import anno
counts_repeats <- cbind(as.data.frame(counts), repeats[, c("gene_id", "transcript_id")])
#count repeat families and classes
counts_sum_repFamily <- aggregate(counts_repeats[,1:(ncol(counts_repeats)-7)], by=list(repFamily=counts_repeats$gene_id), FUN=sum)
counts_sum_repName <- aggregate(cbind(counts_repeats[,1:(ncol(counts_repeats)-7)]), by=list(repFamily=counts_repeats$transcript_id), FUN=sum)
rownames(counts_sum_repFamily) <- counts_sum_repFamily$repFamily
counts_sum_repFamily[,1]<- NULL

#greate dds
dds_repeats <- DESeqDataSetFromMatrix(countData = counts_sum_repFamily, 
                              colData = sample_anno, 
                              design = ~  treatment)

#Estimating size factors
dds_repeats <- estimateSizeFactors(dds_repeats)

#Running the differential expression 
dds_repeats <- DESeq(dds_repeats)
saveRDS(dds_repeats, file.path(results.dir, "dds_repFamily_group.rds"))
