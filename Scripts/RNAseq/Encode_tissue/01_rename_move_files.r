#rename encode files for processing
#load encode metadata_new
meta <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/Encode/data_new/metadata.tsv"))
nrow(meta)
#subset by paired end
meta_sub <- meta[meta$"Run type" =="paired-ended",]
nrow(meta_sub)
#further subset by sequencer
meta_sub <- meta_sub[meta_sub$Platform %in% c("Illumina HiSeq 2000", "Illumina HiSeq 2500"  ),] 
nrow(meta_sub)

#get old names
original_names <- paste0("/omics/groups/OE0219/internal/tinat/Encode/data_new/", meta_sub$"File accession",".fastq.gz")
#check if files exist
table(file.exists(original_names))
#get new names
sra <- meta_sub$dbxrefs
nrow(meta_sub[meta_sub$dbxrefs=="",])
sra <- gsub("SRA:", "", sra)
nrow(meta_sub[meta_sub$dbxrefs=="",])
new_names_temp <- file.path("/omics/groups/OE0219/internal/tinat/Encode/data_new", paste0(sra, "_R", meta_sub$"Paired end",".fastq.gz"))
table(duplicated(new_names_temp))
meta_sub$filename_temp <- new_names_temp
#rename
file.rename(original_names, new_names_temp)
table(file.exists(new_names_temp))

#stratify files by strandedness
meta_sub_split <- split(meta_sub,  meta_sub$"Library strand specific")
length(meta_sub_split)
lapply(meta_sub_split, nrow)

#move files to appropriate folder
library(filesstrings)
for(i in names(meta_sub_split)){
    dir.create(file.path("/omics/groups/OE0219/internal/tinat/Encode/data_new/", i))
    file.move(meta_sub_split[[i]]$filename_temp,file.path("/omics/groups/OE0219/internal/tinat/Encode/data_new", i))
    sra <- meta_sub_split[[i]]$dbxrefs
    sra <- gsub("SRA:", "", sra)
    new_names <- file.path("/omics/groups/OE0219/internal/tinat/Encode/data_new", i,paste0( sra, "_R", meta_sub_split[[i]]$"Paired end",".fastq.gz"))
    meta_sub_split[[i]]$filename <- new_names
}

meta_sub_new <- do.call("rbind", meta_sub_split)
saveRDS(meta_sub_new, file.path("/omics/groups/OE0219/internal/tinat/Encode/data_new","meta_new.rds"))