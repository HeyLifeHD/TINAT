#conda activate R4.0
#R
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#functionx
binnedSum <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
      function(seqname) {
          views <- Views(numvar[[seqname]],
                         bins_per_chrom[[seqname]])
          viewSums(views)
      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

#load data
anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")

#rnaseq -> merged tracks
DMSO_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DMSO.normalized.bigWig")
SB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/SB939.normalized.bigWig")
DAC_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DAC.normalized.bigWig")
DACSB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DACSB939.normalized.bigWig")

#cage
DMSO_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DMSO.bigWig")
SB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/SB939.bigWig")
DAC_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DAC.bigWig")
DACSB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DACSB.bigWig")

track_list <- list(DMSO_rnaseq=DMSO_rnaseq, SB_rnaseq=SB_rnaseq, DAC_rnaseq=DAC_rnaseq, DACSB_rnaseq=DACSB_rnaseq, 
    DMSO_cage=DMSO_cage,SB_cage=SB_cage, DAC_cage=DAC_cage ,DACSB_cage=DACSB_cage)
dir.create("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/")

for(i in names(track_list)){
    temp <- keepStandardChromosomes(track_list[[i]], pruning.mode="coarse")
    mcols(temp)$V1 <- score(temp)
    pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
        paste0("coveragePlot_TSS_",i,".pdf")))
    print(ChIPseeker::plotAvgProf2 (peak = temp ,
                TxDb = anno_gtf, weightCol = "V1"))
    dev.off()
     pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_",i,"_v2.pdf")))
    print(plotPeakProf2(peak = temp, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = anno_gtf, weightCol = "V1",ignore_strand = F))
    dev.off()
    print(i)
}
#combined plot
track_list<- lapply(track_list, function(x){
    x <-  keepStandardChromosomes(x, pruning.mode="coarse")
    mcols(x)$V1 <- score(x)
    x
})

pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","All_combined","_knownGene.pdf")))
plotPeakProf2(peak = track_list, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,
            TxDb = txdb, weightCol = "V1",ignore_strand = F)
dev.off()

#binned sums and then coverage
track_list<- lapply(track_list, function(x){
    #x <- x[seqnames(x)!= "chrM",]
    #x <- dropseqlevels
    x <- keepSeqlevels(x,paste0("chr",c(1:22, "X", "Y")), pruning.mode="coarse")
    x
})
#summarize bigwig
chrlengths <- seqlengths(track_list[[1]])
chrlengths <- GRanges(
  seqnames = names(chrlengths),
 ranges = IRanges(start = rep(1, length(chrlengths)),
                   end = chrlengths
  )
)
tiles <- unlist(tile(range(chrlengths), width=800))
seqlevels(tiles)<- as.character(seqlevels(track_list[[1]]))
track_list_tiled <-   lapply(track_list, function(x){
    score <- coverage(x, weight="score")
    x <- binnedSum(tiles,score,"score")
    x
    }
    )

pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","All_combined","_knownGene_preBinnedGenome.pdf")))
plotPeakProf2(peak = track_list, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,
            TxDb = txdb, weightCol = "V1",ignore_strand = F)
dev.off()

#seperated pltos
#rnaseq
track_list_rnaseq <- list(DMSO_rnaseq=DMSO_rnaseq, SB_rnaseq=SB_rnaseq, DAC_rnaseq=DAC_rnaseq, DACSB_rnaseq=DACSB_rnaseq)
track_list_rnaseq <- lapply(track_list_rnaseq, function(x){
    x <-  keepStandardChromosomes(x, pruning.mode="coarse")
    mcols(x)$V1 <- score(x)
    x
})
pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","RNAseq_combined","_knownGene.pdf")))
print(plotPeakProf2(peak = track_list_rnaseq, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,
            TxDb = txdb, weightCol = "V1",ignore_strand = F))
dev.off()
pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","RNAseq_combined","_knownGene_noScale.pdf")))
print(plotPeakProf2(peak = track_list_rnaseq, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,free_y=FALSE,
            TxDb = txdb, weightCol = "V1",ignore_strand = F))
dev.off()
pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","RNAseq_combined","_deNovo.pdf")))
print(plotPeakProf2(peak = track_list_rnaseq, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,
            TxDb = anno_gtf, weightCol = "V1",ignore_strand = F))
dev.off()
#cage
track_list_cage <- list(DMSO_cage=DMSO_cage,SB_cage=SB_cage, DAC_cage=DAC_cage ,DACSB_cage=DACSB_cage)
track_list_cage <- lapply(track_list_cage, function(x){
    x <-  keepStandardChromosomes(x, pruning.mode="coarse")
    mcols(x)$V1 <- score(x)
    x
})
pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","CAGE_combined","_knownGene.pdf")))
print(plotPeakProf2(peak = track_list_cage, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,
            TxDb = txdb, weightCol = "V1",ignore_strand = F))
dev.off()
pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","CAGE_combined","_knownGene_noScale.pdf")))
print(plotPeakProf2(peak = track_list_cage, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,free_y=FALSE,
            TxDb = txdb, weightCol = "V1",ignore_strand = F))
dev.off()
pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/coverage_plots/",
    paste0("coveragePlot_TSS_","CAGE_combined","_deNovo.pdf")))
print(plotPeakProf2(peak = track_list_cage, upstream = rel(0.2), downstream = rel(0.2),
            conf = 0.95, by = "gene", type = "body", nbin = 800,
            TxDb = anno_gtf, weightCol = "V1",ignore_strand = F))
dev.off()




