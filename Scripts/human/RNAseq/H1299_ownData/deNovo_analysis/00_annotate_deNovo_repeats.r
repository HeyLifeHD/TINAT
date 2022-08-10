#annotate transcripts with repeats
#libraries
library(rtracklayer)
#load data
anno <-  import.gff2("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/mergedTranscripts.gtf")
repeats<- readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))

#distance to closest repeaat based of TSS
dist_N <- distanceToNearest(resize(anno,1), repeats, ignore.strand=TRUE)
anno$nearest_repeat_repName <- NA
anno$nearest_repeat_repClass <- NA
anno$nearest_repeat_repFamily <- NA
anno$dist_nearest_repeat <- NA
anno[queryHits(dist_N),]$nearest_repeat_repName<- as.character(repeats[subjectHits(dist_N), ]$repName)
anno[queryHits(dist_N),]$nearest_repeat_repClass<- as.character(repeats[subjectHits(dist_N), ]$repClass)
anno[queryHits(dist_N),]$nearest_repeat_repFamily<- as.character(repeats[subjectHits(dist_N), ]$repFamily)
anno[queryHits(dist_N),]$dist_nearest_repeat<- as.numeric(dist_N@elementMetadata@listData$distance)

#same for closest LTR12
repeats_sub <- repeats[grep("LTR12", repeats$repName),]
dist_N <- distanceToNearest(resize(anno,1), repeats_sub, ignore.strand=TRUE)
anno$nearest_LTR12repeat_repName <- NA
anno$dist_nearest_LTR12repeat <- NA
anno[queryHits(dist_N),]$nearest_LTR12repeat_repName<- as.character(repeats_sub[subjectHits(dist_N), ]$repName)
anno[queryHits(dist_N),]$dist_nearest_LTR12repeat<- as.numeric(dist_N@elementMetadata@listData$distance)

saveRDS(anno,"/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds" )
