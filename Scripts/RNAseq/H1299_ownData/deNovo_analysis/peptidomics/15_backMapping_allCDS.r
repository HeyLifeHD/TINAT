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

dir.create(output.dir)
#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_gtf <-GenomicFeatures::makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")


#import ORF fasta
#ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs.pep")
#import ORF coordinates
ORF_region_gff <- import.gff3("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs.gff3")
length(ORF_region_gff[ORF_region_gff$type=="CDS"])

#define CDS on genomic coordinates
#get cs regions
cds <- ORF_region_gff[ORF_region_gff$type=="CDS"]
peptide_coord <-list()
#get exon coordinates
exons_transcripts <- GenomicFeatures::exonsBy(anno_gtf, "tx")
names(exons_transcripts)<- GenomicFeatures::transcripts(anno_gtf)$tx_name

#run backmapping
cds_coordinates <- list()
for(i in cds$ID){
    cds_coords <- cds[cds$ID == i,]
    g_coords <- exons_transcripts[[as.character(seqnames(cds[cds$ID == i,]))]]
    cds_coordinates[[i]] <- ensembldb:::.to_genome(g_coords, cds_coords)
    mcols(cds_coordinates[[i]])$accession_new <- cds[cds$ID == i,]$Parent    
    mcols(cds_coordinates[[i]])$transcript_id <- as.character(seqnames(cds[cds$ID == i,]))
    print(round(length(cds_coordinates[[i]]))/(length(cds$ID)), 3)
}
cds_coordinates <- list()

cds_coordinates <- mclapply(cds$ID, function(i){
    cds_coords <- cds[cds$ID == i,]
    g_coords <- exons_transcripts[[as.character(seqnames(cds[cds$ID == i,]))]]
    x <- ensembldb:::.to_genome(g_coords, cds_coords)
    mcols(x)$accession_new <- cds[cds$ID == i,]$Parent    
    mcols(x)$transcript_id <- as.character(seqnames(cds[cds$ID == i,]))
    x
}, mc.cores=10)
saveRDS(cds_coordinates, file.path( "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression", "ORFs_genomic_location.rds"))