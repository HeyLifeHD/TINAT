#locus plots
#libraries
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(Gviz)
#load data
#annotations
#anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
#anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
repeats <- fread("/omics/groups/OE0219/internal/repguide/data/repeats/hg19_repeats.txt")

#rnaseq -> merged tracks
DMSO_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DMSO.bigWig")
SB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/SB939.bigWig")
DAC_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DAC.bigWig")
DACSB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DACSB939.bigWig")

#cage
DMSO_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DMSO.bigWig")
SB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/SB939.bigWig")
DAC_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DAC.bigWig")
DACSB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DACSB.bigWig")

#nanopore
DACSB_nanopore <- import.bw("/omics/groups/OE0219/internal/tinat/210709_nanopore_ashish_stringtie_assembly/pipeline-nanopore-ref-isoforms/bigwig/reads_aln_sorted.bw")

#methylation
DMSO_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
SB_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150389_H2_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
DAC_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
DACSB_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150387_H2_DAC_plus_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")

#H3K9me3
DMSO_9me3 <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150368_DMSO_H3K9me3_ATCACG_L002_R1_complete_filtered.fastq.gz.bigWig")
SB_9me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150384_SB939_H3K9me3_TTAGGC_L002_R1_complete_filtered.fastq.gz.bigWig")
DAC_9me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150352_DAC_H3K9me3_CGATGT_L002_R1_complete_filtered.fastq.gz.bigWig")
DACSB_9me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150360_DACSB_H3K9me3_ACTTGA_L002_R1_complete_filtered.fastq.gz.bigWig")

#H3K4me3
DMSO_4me3 <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150366_DMSO_H3K4me3_ATCACG_L001_R1_complete_filtered.fastq.gz.bigWig")
SB_4me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150382_SB939_H3K4me3_ACTTGA_L001_R1_complete_filtered.fastq.gz.bigWig")
DAC_4me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150350_DAC_H3K4me3_GCCAAT_L001_R1_complete_filtered.fastq.gz.bigWig")
DACSB_4me3 <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150358_DACSB_H3K4me3_CTTGTA_L001_R1_complete_filtered.fastq.gz.bigWig")

#H3K9ac
DMSO_9ac <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150367_DMSO_H3K9ac_CGATGT_L005_R1_complete_filtered.fastq.gz.bigWig")
SB_9ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150383_SB939_H3K9ac_TGACCA_L005_R1_complete_filtered.fastq.gz.bigWig")
DAC_9ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150351_DAC_H3K9ac_TTAGGC_L005_R1_complete_filtered.fastq.gz.bigWig")
DACSB_9ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150359_DACSB_H3K9ac_GATCAG_L005_R1_complete_filtered.fastq.gz.bigWig")

#H3K27ac
DMSO_27ac <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150362_DMSO_H3K27ac_ATCACG_L008_R1_complete_filtered.fastq.gz.bigWig")
SB_27ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150378_SB939_H3K27ac_ACTTGA_L008_R1_complete_filtered.fastq.gz.bigWig")
DAC_27ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150346_DAC_H3K27ac_GCCAAT_L008_R1_complete_filtered.fastq.gz.bigWig")
DACSB_27ac <-import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/GSM2150354_DACSB_H3K27ac_CTTGTA_L008_R1_complete_filtered.fastq.gz.bigWig")

#H3K14ac
chip.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/"
DMSO_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597665_DMSO-7.normalized", ".bigWig")))
SB_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597647_SB-7.normalized", ".bigWig")))
DAC_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597638_DAC-7.normalized", ".bigWig")))
DACSB_14ac <- import.bw(file.path(chip.dir, paste0("GSM2597656_DAC+SB-7.normalized", ".bigWig")))

#H2BK5ac
DMSO_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597670_DMSO-15.normalized", ".bigWig")))
SB_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597652_SB-15.normalized", ".bigWig")))
DAC_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597643_DAC-15.normalized", ".bigWig")))
DACSB_5ac <- import.bw(file.path(chip.dir, paste0("GSM2597661_DAC+SB-15.normalized", ".bigWig")))

#prepare data
#create granges from repeats file
repeats <- makeGRangesFromDataFrame(repeats, seqnames.field="genoName",
    start.field="genoStart", end.field="genoEnd", keep.extra.columns= TRUE)
repeats$repName <- as.character(repeats$repName)

#subset annos for symbol parsing
# anno_original_transcript_symbol <- data.frame(row.names=anno_original[anno_original$type =="transcript",]$transcript_id,
# SYMBOL=  anno_original[anno_original$type =="transcript",]$transcript_name)
anno_transcript_symbol <- data.frame(row.names=anno[anno$type =="transcript",]$transcript_id,
    SYMBOL=  anno[anno$type =="transcript",]$gene_name)
#get class code information
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_transcript_class <- data.frame(row.names=anno_classi$qry_id,
    Class_code=  anno_classi$class_code_simple)
#replace with easy name for gviz function
anno_transcript_class$Class_code <- gsub("non-chimeric (novel)","NonChimeric",anno_transcript_class$Class_code,fixed=TRUE)
anno_transcript_class$Class_code <- gsub("chimeric (novel)","Chimeric",anno_transcript_class$Class_code,fixed=TRUE)
anno_transcript_class$Class_code <- gsub("known","Known",anno_transcript_class$Class_code,fixed=TRUE)

#make granges of methylation
DMSO_meth <- makeGRangesFromDataFrame(DMSO_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(DMSO_meth)<- DMSO_meth$Smoothed_Methylation_Level_H2_DMSO
score(DMSO_meth)<- DMSO_meth$Smoothed_Methylation_Level_H2_DMSO
SB_meth <- makeGRangesFromDataFrame(SB_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(SB_meth)<- SB_meth$Smoothed_Methylation_Level_H2_SB939
score(SB_meth)<- SB_meth$Smoothed_Methylation_Level_H2_SB939
DAC_meth <- makeGRangesFromDataFrame(DAC_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(DAC_meth)<- DAC_meth$Smoothed_Methylation_Level_H2_DAC
score(DAC_meth)<- DAC_meth$Smoothed_Methylation_Level_H2_DAC
DACSB_meth <- makeGRangesFromDataFrame(DACSB_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(DACSB_meth)<- DACSB_meth$Smoothed_Methylation_Level_H2_DAC_plus_SB939 
score(DACSB_meth)<- DACSB_meth$Smoothed_Methylation_Level_H2_DAC_plus_SB939 

#parameters
#create directory
dir.create("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots")
#region
ext <- 5000 
fontSize <- 10
#colors
col <- c("#00AFBB", "#E7B800", "#FC4E07","gray")
names(col) <-c("DAC", "DACandSB939", "DMSO", "SB939")

#select region
temp <- "^DAPK1$"
roi <- anno_original[grep(temp, anno_original$gene_name),]

roi_new <-GRanges(
    seqnames = unique(seqnames(roi)),
    ranges = IRanges(start = min(start(roi)),
                   end =  max(end(roi))))
roi_new$transcript_id <- unique(roi$gene_name)
#plot regions
for (i in 1:length(roi_new)){
    #Define Region
    lim <- c(start(roi_new[i,]), end(roi_new[i,]))
    Chr<- as.character(seqnames(roi_new[i,]))
    range<- GRanges(
    seqnames = Chr,
    ranges = IRanges(start = lim[1]-ext,
                   end = lim[2]+ext))

    #get ideogramm tracks
    itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr,fontcolor="black")

    #genome axis track
    getrack <- GenomeAxisTrack(fontcolor="black")

    #gene annotation
    #original
    # grtrack_original <- GeneRegionTrack(anno_original_gtf, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, 
    #     geneSymbol=TRUE, name="Genecode",collapseTrack=TRUE, fill="darkgray",color="black", cex.feature = 0.5,
    #     fontsize=fontSize-2,fontcolor.title="black", fontcolor.group="black")
    # symbol(grtrack_original) <- as.character(anno_original_transcript_symbol[transcript(grtrack_original),])

    #deNOvo
    grtrack <- GeneRegionTrack(anno_gtf, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, # fill="darkgray",
        geneSymbol=TRUE, name="deNovo",collapseTrack=TRUE,color="black", cex.feature = 0.75,
        fontsize=fontSize,fontcolor.title="black", fontcolor.group="black", Known="darkgray", Chimeric= "orange", NonChimeric="tomato3")
  
    symbol(grtrack) <- as.character(anno_transcript_symbol[transcript(grtrack),])
    feature(grtrack) <- anno_transcript_class[transcript(grtrack),]

    #get annotation tracks of LTR12 repeats
    anno_repeats <- AnnotationTrack(repeats[grep("LTR12", repeats$repName),], name="LTR12", 
        shape = "box",fill="darkgray",color="darkgray", fontsize=fontSize,fontcolor.title="black",fontcolor.feature ="black",
        chromosome=Chr,genome = "hg19", fill="darkgray", featureAnnotation="feature", cex.feature = 0.5)
    feature(anno_repeats)<- repeats[grep("LTR12", repeats$repName),]$repName

    #Data Tracks
    #rnaseq
    #get ylim
    ylim_rnaseq <- c(-(max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,]))))/20), 
        max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,])))))    
    yTicksAt_rnaseq <- c(0, plyr::round_any(max(ylim_rnaseq), accuracy=10, f = floor))
    #prepare tracks
    DMSO_rnaseq_track <- DataTrack(range = DMSO_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black", ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_rnaseq_track <- DataTrack(range = SB_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nSB939", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_rnaseq_track <- DataTrack(range = DAC_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_rnaseq_track <- DataTrack(range = DACSB_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])
   
    #cage
    #get ylim
    ylim_cage <- c(-(max(c(max(score(DMSO_cage[DMSO_cage %over%range,])), max(score(SB_cage[SB_cage %over%range,])), 
        max(score(DAC_cage[DAC_cage %over%range,])), max(score(DACSB_cage[DACSB_cage %over%range,]))))/20), max(c(max(score(DMSO_cage[DMSO_cage %over%range,])), max(score(SB_cage[SB_cage %over%range,])), 
        max(score(DAC_cage[DAC_cage %over%range,])), max(score(DACSB_cage[DACSB_cage %over%range,])))))
    yTicksAt_cage <- c(0, plyr::round_any(max(ylim_cage), accuracy=10, f = floor))
    #prepare tracks
    DMSO_cage_track <- DataTrack(range = DMSO_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_cage_track <- DataTrack(range = SB_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_cage_track <- DataTrack(range = DAC_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage,  yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_cage_track <- DataTrack(range = DACSB_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])
    
    #Methylation
    DMSO_meth_track <- DataTrack(range = DMSO_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05), yTicksAt=c(0,1),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_meth_track <- DataTrack(range = SB_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_meth_track <- DataTrack(range = DAC_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_meth_track <- DataTrack(range = DACSB_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])

    #H3K9me3
    #get ylim
    ylim_9me3 <- c(-(max(c(max(score(DMSO_9me3[DMSO_9me3 %over%range,])), max(score(SB_9me3[SB_9me3 %over%range,])), 
        max(score(DAC_9me3[DAC_9me3 %over%range,])), max(score(DACSB_9me3[DACSB_9me3 %over%range,]))))/20), max(c(max(score(DMSO_9me3[DMSO_9me3 %over%range,])), max(score(SB_9me3[SB_9me3 %over%range,])), 
        max(score(DAC_9me3[DAC_9me3 %over%range,])), max(score(DACSB_9me3[DACSB_9me3 %over%range,])))))
    yTicksAt_9me3 <- c(0, plyr::round_any(max(ylim_9me3), accuracy=10, f = floor))
    #prepare tracks
    DMSO_9me3_track <- DataTrack(range = DMSO_9me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9me3,yTicksAt= yTicksAt_9me3,
        type = c("histogram"), chromosome = Chr, name = "H3K9me3\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_9me3_track <- DataTrack(range = SB_9me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9me3,yTicksAt= yTicksAt_9me3,
        type = c("histogram"), chromosome = Chr, name = "H3K9me3\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_9me3_track <- DataTrack(range = DAC_9me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9me3,yTicksAt= yTicksAt_9me3,
        type = c("histogram"), chromosome = Chr, name = "H3K9me3\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_9me3_track <- DataTrack(range = DACSB_9me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9me3,yTicksAt= yTicksAt_9me3,
        type = c("histogram"), chromosome = Chr, name = "H3K9me3\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])
     
    #H3K4me3
    #get ylim
    ylim_4me3 <- c(-(max(c(max(score(DMSO_4me3[DMSO_4me3 %over%range,])), max(score(SB_4me3[SB_4me3 %over%range,])), 
        max(score(DAC_4me3[DAC_4me3 %over%range,])), max(score(DACSB_4me3[DACSB_4me3 %over%range,]))))/20), max(c(max(score(DMSO_4me3[DMSO_4me3 %over%range,])), max(score(SB_4me3[SB_4me3 %over%range,])), 
        max(score(DAC_4me3[DAC_4me3 %over%range,])), max(score(DACSB_4me3[DACSB_4me3 %over%range,])))))
    yTicksAt_4me3 <- c(0, plyr::round_any(max(ylim_4me3), accuracy=10, f = floor))
    #prepare tracks
    DMSO_4me3_track <- DataTrack(range = DMSO_4me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_4me3, yTicksAt=yTicksAt_4me3, 
        type = c("histogram"), chromosome = Chr, name = "H3K4me3\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_4me3_track <- DataTrack(range = SB_4me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_4me3, yTicksAt=yTicksAt_4me3, 
        type = c("histogram"), chromosome = Chr, name = "H3K4me3\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_4me3_track <- DataTrack(range = DAC_4me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_4me3, yTicksAt=yTicksAt_4me3, 
        type = c("histogram"), chromosome = Chr, name = "H3K4me3\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_4me3_track <- DataTrack(range = DACSB_4me3, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_4me3, yTicksAt=yTicksAt_4me3, 
        type = c("histogram"), chromosome = Chr, name = "H3K4me3\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])

    #H3K9ac
    #get ylim
    ylim_9ac <- c(-(max(c(max(score(DMSO_9ac[DMSO_9ac %over%range,])), max(score(SB_9ac[SB_9ac %over%range,])), 
        max(score(DAC_9ac[DAC_9ac %over%range,])), max(score(DACSB_9ac[DACSB_9ac %over%range,]))))/20), max(c(max(score(DMSO_9ac[DMSO_9ac %over%range,])), max(score(SB_9ac[SB_9ac %over%range,])), 
        max(score(DAC_9ac[DAC_9ac %over%range,])), max(score(DACSB_9ac[DACSB_9ac %over%range,])))))
    yTicksAt_9ac <- c(0, plyr::round_any(max(ylim_9ac), accuracy=10, f = floor))
    #prepare tracks
    DMSO_9ac_track <- DataTrack(range = DMSO_9ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9ac, yTicksAt= yTicksAt_9ac,
        type = c("histogram"), chromosome = Chr, name = "H3K9ac\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_9ac_track <- DataTrack(range = SB_9ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9ac, yTicksAt= yTicksAt_9ac,
        type = c("histogram"), chromosome = Chr, name = "H3K9ac\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_9ac_track <- DataTrack(range = DAC_9ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9ac, yTicksAt= yTicksAt_9ac,
        type = c("histogram"), chromosome = Chr, name = "H3K9ac\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_9ac_track <- DataTrack(range = DACSB_9ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_9ac, yTicksAt= yTicksAt_9ac,
        type = c("histogram"), chromosome = Chr, name = "H3K9ac\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])

    #H3K27ac
    #get ylim
    ylim_27ac <- c(-(max(c(max(score(DMSO_27ac[DMSO_27ac %over%range,])), max(score(SB_27ac[SB_27ac %over%range,])), 
        max(score(DAC_27ac[DAC_27ac %over%range,])), max(score(DACSB_27ac[DACSB_27ac %over%range,]))))/20), max(c(max(score(DMSO_27ac[DMSO_27ac %over%range,])), max(score(SB_27ac[SB_27ac %over%range,])), 
        max(score(DAC_27ac[DAC_27ac %over%range,])), max(score(DACSB_27ac[DACSB_27ac %over%range,])))))
    yTicksAt_27ac <- c(0, plyr::round_any(max(ylim_27ac), accuracy=10, f = floor)/2, plyr::round_any(max(ylim_27ac), accuracy=10, f = floor))
    #prepare tracks
    DMSO_27ac_track <- DataTrack(range = DMSO_27ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_27ac, yTicksAt=yTicksAt_27ac,
        type = c("histogram"), chromosome = Chr, name = "H3K27ac\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_27ac_track <- DataTrack(range = SB_27ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_27ac, yTicksAt=yTicksAt_27ac,
        type = c("histogram"), chromosome = Chr, name = "H3K27ac\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_27ac_track <- DataTrack(range = DAC_27ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_27ac, yTicksAt=yTicksAt_27ac,
        type = c("histogram"), chromosome = Chr, name = "H3K27ac\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_27ac_track <- DataTrack(range = DACSB_27ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_27ac, yTicksAt=yTicksAt_27ac,
        type = c("histogram"), chromosome = Chr, name = "H3K27ac\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])

    #H3K14ac
    #get ylim
    ylim_14ac <- c(-(max(c(max(score(DMSO_14ac[DMSO_14ac %over%range,])), max(score(SB_14ac[SB_14ac %over%range,])), 
        max(score(DAC_14ac[DAC_14ac %over%range,])), max(score(DACSB_14ac[DACSB_14ac %over%range,]))))/20), max(c(max(score(DMSO_14ac[DMSO_14ac %over%range,])), max(score(SB_14ac[SB_14ac %over%range,])), 
        max(score(DAC_14ac[DAC_14ac %over%range,])), max(score(DACSB_14ac[DACSB_14ac %over%range,])))))
    yTicksAt_14ac <- c(0, plyr::round_any(max(ylim_14ac), accuracy=10, f = floor)/2, plyr::round_any(max(ylim_14ac), accuracy=10, f = floor))
    #prepare tracks
    DMSO_14ac_track <- DataTrack(range = DMSO_14ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_14ac, yTicksAt=yTicksAt_14ac,
        type = c("histogram"), chromosome = Chr, name = "H3K14ac\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_14ac_track <- DataTrack(range = SB_14ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_14ac, yTicksAt=yTicksAt_14ac,
        type = c("histogram"), chromosome = Chr, name = "H3K14ac\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_14ac_track <- DataTrack(range = DAC_14ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_14ac, yTicksAt=yTicksAt_14ac,
        type = c("histogram"), chromosome = Chr, name = "H3K14ac\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_14ac_track <- DataTrack(range = DACSB_14ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_14ac, yTicksAt=yTicksAt_14ac,
        type = c("histogram"), chromosome = Chr, name = "H3K14ac\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])

    #H2BK5ac
    #get ylim
    ylim_5ac <- c(-(max(c(max(score(DMSO_5ac[DMSO_5ac %over%range,])), max(score(SB_5ac[SB_5ac %over%range,])), 
        max(score(DAC_5ac[DAC_5ac %over%range,])), max(score(DACSB_5ac[DACSB_5ac %over%range,]))))/20), max(c(max(score(DMSO_5ac[DMSO_5ac %over%range,])), max(score(SB_5ac[SB_5ac %over%range,])), 
        max(score(DAC_5ac[DAC_5ac %over%range,])), max(score(DACSB_5ac[DACSB_5ac %over%range,])))))
    yTicksAt_5ac <- c(0, plyr::round_any(max(ylim_5ac), accuracy=10, f = floor)/2, plyr::round_any(max(ylim_5ac), accuracy=10, f = floor))
    #prepare tracks
    DMSO_5ac_track <- DataTrack(range = DMSO_5ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_5ac, yTicksAt=yTicksAt_5ac,
        type = c("histogram"), chromosome = Chr, name = "H2BK5ac\nDMSO", fill.histogram=col["DMSO"],col.histogram=col["DMSO"])
    SB_5ac_track <- DataTrack(range = SB_5ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_5ac, yTicksAt=yTicksAt_5ac,
        type = c("histogram"), chromosome = Chr, name = "H2BK5ac\nSB", fill.histogram=col["SB939"],col.histogram=col["SB939"])
    DAC_5ac_track <- DataTrack(range = DAC_5ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_5ac, yTicksAt=yTicksAt_5ac,
        type = c("histogram"), chromosome = Chr, name = "H2BK5ac\nDAC", fill.histogram=col["DAC"],col.histogram=col["DAC"])
    DACSB_5ac_track <- DataTrack(range = DACSB_5ac, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_5ac, yTicksAt=yTicksAt_5ac,
        type = c("histogram"), chromosome = Chr, name = "H2BK5ac\nDAC+SB", fill.histogram=col["DACandSB939"],col.histogram=col["DACandSB939"])

     
    #Plot track
    pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots",paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 24)
    plotTracks(list(itrack ,getrack,grtrack, anno_repeats, # grtrack_original,
        DMSO_rnaseq_track,SB_rnaseq_track,DAC_rnaseq_track,DACSB_rnaseq_track,
        DMSO_cage_track,SB_cage_track,DAC_cage_track,DACSB_cage_track,
        DMSO_meth_track,SB_meth_track,DAC_meth_track,DACSB_meth_track,
        DMSO_9me3_track,SB_9me3_track,DAC_9me3_track,DACSB_9me3_track,
        DMSO_4me3_track,SB_4me3_track,DAC_4me3_track,DACSB_4me3_track,
        DMSO_9ac_track,SB_9ac_track,DAC_9ac_track,DACSB_9ac_track,
        DMSO_27ac_track,SB_27ac_track,DAC_27ac_track,DACSB_27ac_track,
        DMSO_14ac_track,SB_14ac_track,DAC_14ac_track,DACSB_14ac_track,
        DMSO_5ac_track,SB_5ac_track,DAC_5ac_track,DACSB_5ac_track
        ),fontcolor.title="black",
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    dev.off()
 
    print(i)
}
