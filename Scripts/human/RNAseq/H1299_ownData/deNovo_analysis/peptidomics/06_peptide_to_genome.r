#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)
library(Biostrings)
library(rtracklayer)
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"
dir.create(output.dir)
#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
repeats<- readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))

#get transcripts of anno
anno_transcripts <- 
#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
comp_col <- col[c(2,8,6)]
names(comp_col)<-c( "DAC_vs_DMSO", "SB939_vs_DMSO",  "DACandSB939_vs_DMSO")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
ltr_col <- brewer.pal(5, "Accent")
names(ltr_col)<- c("LTR12","LTR12C","LTR12D","LTR12F", "no LTR12")
#peptidomics list
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))

#import ORF fasta
ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")
#import ORF coordinates
ORF_region_gff <- import.gff3("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs.gff3")

#subset peptides that originatee from our orf list
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
#select 2/3 candidates
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]$sequences))
peptides_new_ORF_oi <- peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]
peptides_new_ORF_oi[peptides_new_ORF_oi$accession_new ==peptides_new_ORF_oi$accession_new[duplicated(peptides_new_ORF_oi$accession_new)],]  
 
$accession_new)   ]        
#check for peptide in the respective orf
#renmae orfs
ORF_names_old <- names(ORFs)
names(ORF_names_old)<-sapply(strsplit(names(ORFs), " type"), "[",1)
names(ORFs)<- sapply(strsplit(names(ORFs), " type"), "[",1)
table(peptides_new_ORF_oi$accession_new %in% names(ORFs))
#make granges of ORF position
temp <- sapply(strsplit(ORF_names_old, ":",fixed=TRUE), "[",5)
temp <- sapply(strsplit(temp, ")",fixed=TRUE), "[",1)
strand <- sapply(strsplit(temp, "(",fixed=TRUE), "[",2)
range <- sapply(strsplit(temp, "(",fixed=TRUE), "[",1)
ORF_regions <-IRanges(start =as.numeric(sapply(strsplit(range, "-",fixed=TRUE), "[",1)),
                    end =  as.numeric(sapply(strsplit(range, "-",fixed=TRUE), "[",2)))
names(ORF_regions)<- names(ORFs)   

#loop over peptides_new_ORF_oi and find pattern match of ORFs
peptide_regions <- ORF_regions[peptides_new_ORF_oi$accession_new,]
pattern_match <- list()
peptide_regions_new <- IRangesList()
for(i in 1: length(peptides_new_ORF_oi$accession_new)){
    accession_oi <- peptides_new_ORF_oi[i,]$accession_new
    sequence_oi <- peptides_new_ORF_oi[i,]$sequence 
    pattern_match[[i]] <- vmatchPattern(sequence_oi, ORFs[accession_oi,])
    #get pattern in relationship to ORF
    peptide_regions_new[[i]] <- IRanges(start = as.numeric(start(peptide_regions[accession_oi,])+((start(pattern_match[[i]])-1)*3)),
                    end =  as.numeric(start(peptide_regions[accession_oi,])+((start(pattern_match[[i]])-1)*3)+((width(peptides_new_ORF_oi[i,]$sequences)*3)-1)))
                    #end =  as.numeric(start(peptide_regions[accession_oi,])+((start(pattern_match[[i]])-1)*3)+((8*3)-1)))
} 
table(sapply(peptide_regions_new, width)/3 ==width(peptides_new_ORF_oi$sequences))
table(length(peptide_regions_new)==length(peptides_new_ORF_oi$accession_new))

peptide_regions_new <- unlist(peptide_regions_new)
mcols(peptide_regions_new)$accession_new <-peptides_new_ORF_oi$accession_new 
mcols(peptide_regions_new)$sequences<-peptides_new_ORF_oi$sequences 
temp <-sapply(strsplit(mcols(peptide_regions_new)$accession_new , ".p"), "[",1)
mcols(peptide_regions_new)$transcript_id<- gsub("ORF_", "",temp)

#get exon coordinates
exons_transcripts <- exonsBy(anno_gtf, "tx")
names(exons_transcripts)<- transcripts(anno_gtf)$tx_name
exons_transcripts_sub <- exons_transcripts[mcols(peptide_regions_new)$transcript_id]

#define regions for an example
peptide_coord <-list()
for(i in 1:length(peptide_regions_new)){
    cds_coords <- peptide_regions_new[i,]
    g_coords <- exons_transcripts_sub[[mcols(peptide_regions_new[i,])$transcript_id]]
    peptide_coord[[i]] <- ensembldb:::.to_genome(g_coords, cds_coords)
    mcols(peptide_coord[[i]])$accession_new <- mcols(peptide_regions_new[i,])$accession_new
    mcols(peptide_coord[[i]])$transcript_id <- mcols(peptide_regions_new[i,])$transcript_id
    mcols(peptide_coord[[i]])$sequences <- mcols(peptide_regions_new[i,])$sequences
}

#annotate granges list top ltr sequences
peptide_coord <- lapply(peptide_coord, function(x){
    dist_N <- distanceToNearest(x, repeats_sub, ignore.strand=FALSE)
    mcols(x)$nearest_LTR12repeat_repName <- NA
    mcols(x)$dist_nearest_LTR12repeat <- NA
    mcols(x[queryHits(dist_N),])$nearest_LTR12repeat_repName<- as.character(repeats_sub[subjectHits(dist_N), ]$repName)
    mcols(x[queryHits(dist_N),])$dist_nearest_LTR12repeat<- as.numeric(dist_N@elementMetadata@listData$distance)
    mcols(x)$dist_nearest_LTR12repeat <- ifelse(is.na(mcols(x)$dist_nearest_LTR12repeat), Inf,mcols(x)$dist_nearest_LTR12repeat )
    x
})

#unlist for orf export
peptide_coord_unlist <- unlist(GRangesList(peptide_coord))
export.bed(peptide_coord_unlist,file.path(output.dir, "2of3DACSB_unique","peptide_coordinates.bed") )
saveRDS(peptide_coord,file.path(output.dir,"2of3DACSB_unique", "peptide_coordinates.rds") )

#calculate number of ltr12 ovelrapping peptides
ltr12_overlap <- as.data.frame(table(sapply(peptide_coord, function(x){
    x <- x$dist_nearest_LTR12repeat ==0
ifelse(length(x)>1, ifelse(sum(x)==1,TRUE,FALSE), x)
})))

ltr12_type <- as.data.frame(table(sapply(peptide_coord, function(x){
    x <-ifelse(x$dist_nearest_LTR12repeat ==0, x$nearest_LTR12repeat_repName, "no LTR12")
    ifelse(length(x)>1, unique(x),x)
})))

labs <- paste0(ltr12_type$Freq)
pie<- ggpie(ltr12_type, "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=ltr_col,
                fill = "Var1", color = "black")+rremove("legend.title")
pdf(file.path(output.dir, "2of3DACSB_unique", "peptide_backmap_LTR12anno.pdf"), height=5, width=5)
pie
dev.off()