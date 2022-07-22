#plot feature counts of repeats -SRP230446 - paired end
#libraries
library(rtracklayer)
library(data.table)
library(DESeq2)

#folders
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_FeatureCountsRepeats_analysis/"
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
counts <- fread("/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_FeatureCountsRepeats/repeats_featureCounts.txt")
anno <- counts[,1:6]
counts <- counts[,7:ncol(counts)]
#rename counts sample names
#create annotation sheet

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
temp <- sample_anno
temp$Sample_ID <- gsub(".", "-", sample_anno$Sample_ID, fixed=TRUE)
temp <- temp[order(temp$Sample_ID),]
write.table(temp, file.path(results.dir, "sample_anno.csv"), sep=";", row.names=FALSE, col.names=TRUE, quote=FALSE)
rownames(sample_anno)<-  as.character(sample_anno$Sample_ID)
#rename counts
colnames(counts) <- sapply(strsplit(colnames(counts),"/"),"[",11)
colnames(counts) <-sapply(strsplit(colnames(counts),".sorted"),"[",1)
colnames(counts) <- gsub("-",".", colnames(counts), fixed=TRUE)
rownames(sample_anno) %in% colnames(counts)

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
                              design =  ~ cellline + treatment )

#Estimating size factors
dds_repeats <- estimateSizeFactors(dds_repeats)

#Running the differential expression 
dds_repeats <- DESeq(dds_repeats)
saveRDS(dds_repeats, file.path(results.dir, "dds_repFamily_group.rds"))
