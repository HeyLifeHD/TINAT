#locus plots
#libraries
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(Gviz)
library(BSgenome.Mmusculus.UCSC.mm10)

#load data
#annotations
#anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.gtf")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.mergedTranscripts.gtf.tmap"))
anno_original <-  import.gff2("/omics/groups/OE0219/internal/genomes/Mmusculus/mm10/rnaseq/gencode.vM19.annotation.gtf")
repeats <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/","data","repeats_mm10.rds"))

#rnaseq -> merged tracks treatment experiment
DMSO_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/DMSO.normalized.bigWig")
SB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/SB939.normalized.bigWig")
DAC_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/DAC.normalized.bigWig")
DACSB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/DACSB939.normalized.bigWig")

#rnaseq -> merged tracks setb1kd experiment
SETB1KD_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/bigwig/SETDB1_KO.bigWig")
CONTROL_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/bigwig/control.bigWig")

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_analysis"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Read in Data
DEG_results_list_DACSB<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#setb1
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220811_RNAseqSETB1KD_deNovoB16DACSB_analysis"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Read in Data
DEG_results_list_SETB1<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#prepare data
#create granges from repeats file
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

#parameters
#region
ext <- 1000
fontSize <- 10
#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
cond_col <- col[c(4,6)]
names(cond_col)<-c("control",  "SETB1_KD")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")



#with nanopore data or alignment data
dir.create("/omics/groups/OE0219/internal/tinat/mouse_project/locus_plot")
dir.create("/omics/groups/OE0219/internal/tinat/mouse_project/locus_plot/B16_treat_KD")
dir.create("/omics/groups/OE0219/internal/tinat/mouse_project/locus_plot/B16_treat_KD_withAlign")

#select region
#subset upregulated genes
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff
DEG_results_list_DACSB_up_id <- lapply(DEG_results_list_DACSB, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange>lfc),]
    x$transcript_id
})
DEG_results_list_SETB1_up_id <- lapply(DEG_results_list_SETB1, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange>lfc),]
    x$transcript_id
})
transcript_oi <- DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO[DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO %in%  DEG_results_list_SETB1_up_id$SETDB1_vs_control]
anno_transcript <- anno[anno$type == "transcript",]
roi_new <- anno_transcript
roi_new <- roi_new[roi_new$transcript_id %in% transcript_oi,]

#roi based on coordinates 
# roi_new <- GRanges(
#     seqnames = c("chr8", "chr8", "chr16", "chr10", "chr20","chr22","chr22", "chr9", "chr9","chr12", "chr1", "chr8", "chr19" ),
#     ranges = IRanges(start = c(101267295, 101267295,20942433,114128500,55922500,30991718,30991718,91973702, 91973702, 7144000,119713854, 77347011,8443894),
#                    end =c(101350432,101330400, 21172762,114190134,55956258,31039779,31027974,92112571,92116000,7179000,119727577,77361001,8454545)))
# roi_new$transcript_id <- c("RNF19A_1","RNF19A_2", "DNAH3", "ACSL5" , "RAE1","TCN2", "TCN2","SEMA4D_1","SEMA4D_2" ,"C1S"  ,"MSTRG.1804","MSTRG.32198","MSTRG.16351"   )
# roi_new$transcript_id <- paste0(roi_new$transcript_id ,"_selected")

#plot regions
for (i in 1:length(roi_new)){
    tryCatch({
    #Define Region
    lim <- c(start(roi_new[i,]), end(roi_new[i,]))
    Chr<- as.character(seqnames(roi_new[i,]))
    range<- GRanges(
    seqnames = Chr,
    ranges = IRanges(start = lim[1]-ext,
                   end = lim[2]+ext))

    #get ideogramm tracks
    itrack <- IdeogramTrack(genome = "mm10", chromosome =Chr,fontcolor="black")

    #genome axis track
    getrack <- GenomeAxisTrack(fontcolor="black")

    #deNOvo anno
    grtrack <- GeneRegionTrack( , chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, # fill="darkgray",
        geneSymbol=TRUE, name="deNovo",collapseTrack=TRUE,color="black", cex.feature = 0.75,
        fontsize=fontSize,fontcolor.title="black", fontcolor.group="black", 
        Known=class_col["known"], Chimeric=class_col["chimeric (novel)"], NonChimeric=class_col["non-chimeric (novel)"])
    symbol(grtrack) <- as.character(anno_transcript_symbol[transcript(grtrack),])
    feature(grtrack) <- anno_transcript_class[transcript(grtrack),]

    #get annotation tracks of repeats
    anno_repeats <- AnnotationTrack(repeats, name="repeats", 
        shape = "box",fill="darkgray",color="darkgray", fontsize=fontSize,fontcolor.title="black",fontcolor.feature ="black",
        chromosome=Chr,genome = "mm10", fill="darkgray", featureAnnotation="feature", cex.feature = 0.5)
    feature(anno_repeats)<- repeats$repName

    #Data Tracks
    #rnaseq treatment
    #get ylim
    ylim_rnaseq <- c(-(max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,]))))/20), 
        max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,])))))    
    yTicksAt_rnaseq <- c(0, plyr::round_any(max(ylim_rnaseq), accuracy=10, f = floor))

    #prepare tracks
    DMSO_rnaseq_track <- DataTrack(range = DMSO_rnaseq, genome = "mm10",
        fontsize=fontSize,fontcolor.title="black",col.axis="black", ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDMSO", fill.histogram=treat_col["DMSO"],col.histogram=treat_col["DMSO"])
    SB_rnaseq_track <- DataTrack(range = SB_rnaseq, genome = "mm10",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nSB939", fill.histogram=treat_col["SB939"],col.histogram=treat_col["SB939"])
    DAC_rnaseq_track <- DataTrack(range = DAC_rnaseq, genome = "mm10",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC", fill.histogram=treat_col["DAC"],col.histogram=treat_col["DAC"])
    DACSB_rnaseq_track <- DataTrack(range = DACSB_rnaseq, genome = "mm10",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])

    #rnaseq knockdown
    #get ylim
    ylim_rnaseq_cond <- c(-(max(c(max(score(SETB1KD_rnaseq[SETB1KD_rnaseq %over%range,])), max(score(CONTROL_rnaseq[CONTROL_rnaseq %over%range,]))
       ))/20), 
        max(c(max(score(SETB1KD_rnaseq[SETB1KD_rnaseq %over%range,])), max(score(CONTROL_rnaseq[CONTROL_rnaseq %over%range,]))
        )))    
    yTicksAt_rnaseq_cond <- c(0, plyr::round_any(max(ylim_rnaseq_cond), accuracy=10, f = floor))
    #prepare tracks
    SETB1KD_rnaseq_track <- DataTrack(range = SETB1KD_rnaseq, genome = "mm10",
        fontsize=fontSize,fontcolor.title="black",col.axis="black", ylim=ylim_rnaseq_cond, yTicksAt=yTicksAt_rnaseq_cond,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nSETB1-KD", fill.histogram=cond_col["SETB1_KD"],col.histogram=treat_col["SETB1_KD"])
    CONTROL_rnaseq_track <- DataTrack(range = CONTROL_rnaseq, genome = "mm10",
        fontsize=fontSize,fontcolor.title="black",col.axis="black", ylim=ylim_rnaseq_cond, yTicksAt=yTicksAt_rnaseq_cond,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\ncontrol", fill.histogram=cond_col["control"],col.histogram=treat_col["control"])

    #bam track for dacsb
    DAC_SB_alignment <- AlignmentsTrack("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/HISAT2/aligned_sorted/AS-811016-LR-62189.sorted.bam",
        isPaired = TRUE, type=c("sashimi"),col.sashimi=treat_col["DACandSB939"],sashimiHeight=.5,max.height=.5,
        fontsize=fontSize,fontcolor.title="black",col.axis="black")#,  sashimiFilter =  intronicParts(anno_gtf) ,sashimiFilterTolerance = 10000L
    
    #Plot track
    pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/locus_plot/B16_treat_KD",
        paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 8)
    print(plotTracks(list(itrack ,getrack,anno_repeats,grtrack,# grtrack_original,
        DMSO_rnaseq_track,SB_rnaseq_track,DAC_rnaseq_track,DACSB_rnaseq_track,
        CONTROL_rnaseq_track, SETB1KD_rnaseq_track
        ),fontcolor.title="black",sizes=c(0.5,1,0.5,2,  1,1,1,1, 1,1),
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr))
    dev.off()

    pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/locus_plot/B16_treat_KD_withAlign",
        paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 8)
    print(plotTracks(list(itrack ,getrack,anno_repeats,grtrack, # grtrack_original,
        DMSO_rnaseq_track,SB_rnaseq_track,DAC_rnaseq_track,DACSB_rnaseq_track,DAC_SB_alignment,
        CONTROL_rnaseq_track, SETB1KD_rnaseq_track
        ),fontcolor.title="black",sizes=c(0.5,1,0.5,2,  1,1,1,1,  2.5,  1,1  ),
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr))
    dev.off()
    print(i)

    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

