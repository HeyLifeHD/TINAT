#load encode metadata
meta <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/Encode/data/metadata.csv"))
meta_sub <- meta[meta$"Run type" =="paired-ended",]
meta_split <- split(meta_sub, meta_sub$dbxrefs)

meta_new <- list()
for(i in unique(meta_sub$dbxrefs)[unique(meta_sub$dbxrefs)!=""]){
    meta_new[[i]] <- data.frame(sample=i, fastq_1=paste0("/omics/groups/OE0219/internal/tinat/Encode/data/",meta_split[[i]][meta_split[[i]]$"Paired end"==1,]$"File accession", ".fastq.gz"),
    fastq_2=paste0("/omics/groups/OE0219/internal/tinat/Encode/data/",meta_split[[i]][meta_split[[i]]$"Paired end"==2,]$"File accession", ".fastq.gz"),
    strandedness=unique(meta_split[[i]]$"Library strand specific"))
}
meta_new <- do.call("rbind",meta_new)
dir.create("/omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing")
write.table(meta_new,"/omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing/design.csv" , sep=",",  row.names=F, col.names=TRUE)