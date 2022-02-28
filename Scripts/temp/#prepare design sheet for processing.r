#prepare design sheet for processing
library(data.table)
meta <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/fastq/meta_data.txt")

#subset chip
meta <- meta[meta$LibrarySelection =="ChIP",]
meta_sub <- meta[meta$"DATASTORE filetype" =="fastq,sra",]


unique(meta[meta$"DATASTORE filetype" =="fastq,sra",]$Antibody)
 unique(meta[meta$"DATASTORE filetype" =="bam,sra",]$Antibody)
paste0(meta_sub$treatment, "_", meta_sub$Antibody)
group,replicate,fastq_1,fastq_2,antibody,control