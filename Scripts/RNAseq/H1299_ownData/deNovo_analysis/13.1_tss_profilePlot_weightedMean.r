#heatmap of TSS
#libraries
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(peakSeason)

#include weighted mean as possibility to normalize data in plot

plot_profile2 = function(mat_list, summarizeBy = 'mean', color = NULL, ci = FALSE,weightedMean =FALSE,
                        legend_font_size = 1, condition = NULL, collapse_reps = TRUE, line_size = 2, axis_lwd = 3, axis_font_size = 1,
                        sample_names = NULL, auto_y_axis = TRUE, y_lims = NULL){

  coldata = mat_list$cdata

  if(!is.null(condition)){
    if(!condition %in% colnames(coldata)){
      warning(paste0(condition, " does not exists in coldata. Here are available columns."))
      print(coldata)
      stop()
    }else{
      colnames(coldata)[which(colnames(coldata) == condition)] = "group_condition"
      coldata$group_condition = as.character(coldata$group_condition)
    }
    sample_conditons = as.character(coldata$group_condition)
  }else{
    sample_conditons = NULL
  }

  if(ci){
    ci_list =  peakSeason:::estimateCI(mats = mat_list, group = sample_conditons, collapse_reps = collapse_reps)
  }

  mat_list = peakSeason:::summarizeMats(mats = mat_list, summarizeBy = summarizeBy, group = sample_conditons, collapse_reps = collapse_reps)



  size = as.character(mat_list$param["size"])
  mat_list = mat_list[1:(length(mat_list)-1)]

  #normalize by weighted mean 
    if(weightedMean){
      mat_list =lapply(mat_list, function(x){
          x <- x/weighted.mean(x)
      })
  }

  if(is.null(condition)){
    names(mat_list) = coldata$sample_names
  }

  #return(mat_list)

  if(is.null(color)){
    color = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:length(mat_list)]
    #names(color) = names(mat_list)
  }else{
    if(length(color) != length(mat_list)){
      warning("Insufficient number of colors. Randomly choosing few colors")
      color = c(color,
                sample(colors(), size = length(mat_list) - length(color), replace = FALSE))
    }
  }

  if(is.null(condition)){
    names(color) = names(mat_list)
  }else{
    if(collapse_reps){
      names(color) = names(mat_list)
    }else{
      cond = unique(sample_conditons)
      cond_color = color[1:length(cond)]
      names(cond_color) = cond
      print(cond_color)
      color = cond_color[sample_conditons]
    }
  }


  if(auto_y_axis){
    yl = round(c(min(unlist(lapply(mat_list, min, na.rm = TRUE))),
           max(unlist(lapply(mat_list, max, na.rm = TRUE)))), digits = 2)
  }else{
    yl = round(c(0, max(unlist(lapply(mat_list, max, na.rm = TRUE)))), digits = 2)
  }

  if(ci){
    yl[length(yl)] = yl[length(yl)] + max(unlist(lapply(ci_list, max, na.rm = TRUE)))
  }

  if(auto_y_axis){
    yl[1] = round(min(yl) - (min(yl)*0.10), digits = 2)
    yl[length(yl)] = round(max(yl) + (max(yl)*0.10), digits = 2)
  }

  if(!is.null(y_lims)){
    if(length(y_lims) < 2){
      warning("y_lim requires lower and upper limits. Auto adjusting y axis limits")
      yl[1] = round(min(yl) - (min(yl)*0.10), digits = 2)
      yl[length(yl)] = round(max(yl) + (max(yl)*0.10), digits = 2)
    }else{
      yl = sort(as.numeric(as.character(y_lims)))[1:2]
    }
  }

  xlabs = c(sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 1), 0,
            sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 2))
  xticks = c(0,
             as.integer(length(mat_list[[1]])/sum(as.numeric(xlabs[1]), as.numeric(xlabs[3])) * as.numeric(xlabs[1])),
             length(mat_list[[1]]))

  # mat_list = data.table::rbindlist(lapply(mat_list, function(x){
  #   x = data.table::melt(data = x)
  #   data.table::setDT(x = x)
  #   x[,xax := 1:nrow(x)]
  #   x
  # }), idcol = "sample")


  # g = ggplot(data = mat_list, aes(x = xax, y = value, col = sample))+geom_line()+
  #   scale_y_continuous(breaks = yl, labels = yl, expand = c(0, 0), limits = yl)+
  #   scale_x_continuous(breaks = xticks, labels = c(paste0("-", xlabs[1]), xlabs[2], xlabs[3]),
  #                      expand = c(0, 0), limits = xticks[c(1, 3)])+
  #   theme(panel.background  = element_blank(), axis.line.x = element_line(color="black", size = 0.2),
  #         axis.line.y = element_line(color="black", size = 0.2))
  #
  # return(g)


  par(mar = c(3, 3, 4, 2), font = 4)
  plot(mat_list[[1]], frame.plot = FALSE, axes = F,
       xlab = '', ylab = '', type = 'l',
       lwd = line_size,
       col = color[1],#color[names(mat_list)[[1]]]
       ylim = c(min(yl), max(yl)))

  if(ci){
    polygon(x = c(1:length(mat_list[[1]]), rev(1:length(mat_list[[1]]))),
            y = c(mat_list[[1]]-ci_list[[1]], rev(mat_list[[1]]+ci_list[[1]])),
            col = grDevices::adjustcolor(col = color[1], #color[names(mat_list)[[1]]],
                                         alpha.f = 0.4), border = NA)
  }
  if(length(mat_list) > 1){
    for(i in 2:length(mat_list)){

      points(mat_list[[i]], type = 'l', lwd = line_size, col = color[i])#color[names(mat_list)[[i]]])
      if(ci){
        polygon(x = c(1:length(mat_list[[i]]), rev(1:length(mat_list[[i]]))),
                y = c(mat_list[[i]]-ci_list[[i]], rev(mat_list[[i]]+ci_list[[i]])),
                col = grDevices::adjustcolor(col = color[i], #color[names(mat_list)[[i]]],
                                             alpha.f = 0.4), border = NA)
      }
    }
  }

  xlabs = c(sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 1), 0,
            sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 2))
  xticks = c(0,
             as.integer(length(mat_list[[1]])/sum(as.numeric(xlabs[1]), as.numeric(xlabs[3])) * as.numeric(xlabs[1])),
             length(mat_list[[1]]))

  axis(side = 1, at = xticks,
       labels = c(paste0("-", xlabs[1]), xlabs[2], xlabs[3]), lty = 1, lwd = axis_lwd,
       font = 2, cex.axis = axis_font_size, line = 0.5)
  axis(side = 2, at = yl, labels = yl, lty = 1, lwd = axis_lwd, las = 2, font = 2, cex = axis_font_size)

  if(is.null(condition)){
    peakSeason:::add_legend(x = "topright", bty = "n", legend = names(color),
               col = color, lty = 1, lwd = 3, cex = legend_font_size, ncol = 2)
  }else{
    if(collapse_reps){
      peakSeason:::add_legend(x = "topright", bty = "n", legend = names(color),
                 col = color, lty = 1, lwd = 3, cex = legend_font_size, ncol = 2)
    }else{
      peakSeason:::add_legend(x = "topright", bty = "n", legend = names(cond_color),
                 col = cond_color, lty = 1, lwd = 3, cex = legend_font_size, ncol = 2)
    }
  }
}

plot_profile3 = function(mat_list, summarizeBy = 'mean', color = NULL, ci = FALSE,weightedMean =FALSE,
                        legend_font_size = 1, condition = NULL, collapse_reps = TRUE, line_size = 2, axis_lwd = 3, axis_font_size = 1,
                        sample_names = NULL, auto_y_axis = TRUE, y_lims = NULL){

  coldata = mat_list$cdata

  if(!is.null(condition)){
    if(!condition %in% colnames(coldata)){
      warning(paste0(condition, " does not exists in coldata. Here are available columns."))
      print(coldata)
      stop()
    }else{
      colnames(coldata)[which(colnames(coldata) == condition)] = "group_condition"
      coldata$group_condition = as.character(coldata$group_condition)
    }
    sample_conditons = as.character(coldata$group_condition)
  }else{
    sample_conditons = NULL
  }

  if(ci){
    ci_list =  peakSeason:::estimateCI(mats = mat_list, group = sample_conditons, collapse_reps = collapse_reps)
  }

  mat_list = peakSeason:::summarizeMats(mats = mat_list, summarizeBy = summarizeBy, group = sample_conditons, collapse_reps = collapse_reps)



  size = as.character(mat_list$param["size"])
  mat_list = mat_list[1:(length(mat_list)-1)]

  #normalize by weighted mean 
    if(weightedMean){
      mat_list =lapply(mat_list, function(x){
          x[is.na(x)]<-0
          x <- x/weighted.mean(x)
          x
      })
  }

  if(is.null(condition)){
    names(mat_list) = coldata$sample_names
  }

  #return(mat_list)

  if(is.null(color)){
    color = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:length(mat_list)]
    #names(color) = names(mat_list)
  }else{
    if(length(color) != length(mat_list)){
      warning("Insufficient number of colors. Randomly choosing few colors")
      color = c(color,
                sample(colors(), size = length(mat_list) - length(color), replace = FALSE))
    }
  }

  if(is.null(condition)){
    names(color) = names(mat_list)
  }else{
    if(collapse_reps){
      names(color) = names(mat_list)
    }else{
      cond = unique(sample_conditons)
      cond_color = color[1:length(cond)]
      names(cond_color) = cond
      print(cond_color)
      color = cond_color[sample_conditons]
    }
  }


  if(auto_y_axis){
    yl = round(c(min(unlist(lapply(mat_list, min, na.rm = TRUE))),
           max(unlist(lapply(mat_list, max, na.rm = TRUE)))), digits = 2)
  }else{
    yl = round(c(0, max(unlist(lapply(mat_list, max, na.rm = TRUE)))), digits = 2)
  }

  if(ci){
    yl[length(yl)] = yl[length(yl)] + max(unlist(lapply(ci_list, max, na.rm = TRUE)))
  }

  if(auto_y_axis){
    yl[1] = round(min(yl) - (min(yl)*0.10), digits = 2)
    yl[length(yl)] = round(max(yl) + (max(yl)*0.10), digits = 2)
  }

  if(!is.null(y_lims)){
    if(length(y_lims) < 2){
      warning("y_lim requires lower and upper limits. Auto adjusting y axis limits")
      yl[1] = round(min(yl) - (min(yl)*0.10), digits = 2)
      yl[length(yl)] = round(max(yl) + (max(yl)*0.10), digits = 2)
    }else{
      yl = sort(as.numeric(as.character(y_lims)))[1:2]
    }
  }

  xlabs = c(sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 1), 0,
            sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 2))
  xticks = c(0,
             as.integer(length(mat_list[[1]])/sum(as.numeric(xlabs[1]), as.numeric(xlabs[3])) * as.numeric(xlabs[1])),
             length(mat_list[[1]]))

  # mat_list = data.table::rbindlist(lapply(mat_list, function(x){
  #   x = data.table::melt(data = x)
  #   data.table::setDT(x = x)
  #   x[,xax := 1:nrow(x)]
  #   x
  # }), idcol = "sample")


  # g = ggplot(data = mat_list, aes(x = xax, y = value, col = sample))+geom_line()+
  #   scale_y_continuous(breaks = yl, labels = yl, expand = c(0, 0), limits = yl)+
  #   scale_x_continuous(breaks = xticks, labels = c(paste0("-", xlabs[1]), xlabs[2], xlabs[3]),
  #                      expand = c(0, 0), limits = xticks[c(1, 3)])+
  #   theme(panel.background  = element_blank(), axis.line.x = element_line(color="black", size = 0.2),
  #         axis.line.y = element_line(color="black", size = 0.2))
  #
  # return(g)


  par(mar = c(3, 3, 4, 2), font = 4)
  plot(mat_list[[1]], frame.plot = FALSE, axes = F,
       xlab = '', ylab = '', type = 'l',
       lwd = line_size,
       col = color[1],#color[names(mat_list)[[1]]]
       ylim = c(min(yl), max(yl)))

  if(ci){
    polygon(x = c(1:length(mat_list[[1]]), rev(1:length(mat_list[[1]]))),
            y = c(mat_list[[1]]-ci_list[[1]], rev(mat_list[[1]]+ci_list[[1]])),
            col = grDevices::adjustcolor(col = color[1], #color[names(mat_list)[[1]]],
                                         alpha.f = 0.4), border = NA)
  }
  if(length(mat_list) > 1){
    for(i in 2:length(mat_list)){

      points(mat_list[[i]], type = 'l', lwd = line_size, col = color[i])#color[names(mat_list)[[i]]])
      if(ci){
        polygon(x = c(1:length(mat_list[[i]]), rev(1:length(mat_list[[i]]))),
                y = c(mat_list[[i]]-ci_list[[i]], rev(mat_list[[i]]+ci_list[[i]])),
                col = grDevices::adjustcolor(col = color[i], #color[names(mat_list)[[i]]],
                                             alpha.f = 0.4), border = NA)
      }
    }
  }

  xlabs = c(sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 1), 0,
            sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 2))
  xticks = c(0,
             as.integer(length(mat_list[[1]])/sum(as.numeric(xlabs[1]), as.numeric(xlabs[3])) * as.numeric(xlabs[1])),
             length(mat_list[[1]]))

  axis(side = 1, at = xticks,
       labels = c(paste0("-", xlabs[1]), xlabs[2], xlabs[3]), lty = 1, lwd = axis_lwd,
       font = 2, cex.axis = axis_font_size, line = 0.5)
  axis(side = 2, at = yl, labels = yl, lty = 1, lwd = axis_lwd, las = 2, font = 2, cex = axis_font_size)

  if(is.null(condition)){
    peakSeason:::add_legend(x = "topright", bty = "n", legend = names(color),
               col = color, lty = 1, lwd = 3, cex = legend_font_size, ncol = 2)
  }else{
    if(collapse_reps){
      peakSeason:::add_legend(x = "topright", bty = "n", legend = names(color),
                 col = color, lty = 1, lwd = 3, cex = legend_font_size, ncol = 2)
    }else{
      peakSeason:::add_legend(x = "topright", bty = "n", legend = names(cond_color),
                 col = cond_color, lty = 1, lwd = 3, cex = legend_font_size, ncol = 2)
    }
  }
}

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"
input.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/6_tracks_normalized"
#load data
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#Set specifications 
cutoff <- 0.01 #set FDR cutoff
l2fc <- 2##set logfold2 cutoff

#annotations
#anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
#anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
repeats <- fread("/omics/groups/OE0219/internal/repguide/data/repeats/hg19_repeats.txt")
responsive_TSS <- readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/responsive_TSSs.RDS")

#rnaseq
bw_rna <- list("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DMSO.normalized.bigWig","/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/SB939.normalized.bigWig",
    "/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DAC.normalized.bigWig","/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DACSB939.normalized.bigWig")

#cage
bw_cage <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DMSO.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/SB939.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DAC.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DACSB.bigWig")

#H3K9me3
bw_9me3 <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150368_DMSO_H3K9me3_ATCACG_L002_R1_complete_filtered.fastq.gz.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150384_SB939_H3K9me3_TTAGGC_L002_R1_complete_filtered.fastq.gz.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150352_DAC_H3K9me3_CGATGT_L002_R1_complete_filtered.fastq.gz.bigWig","/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150360_DACSB_H3K9me3_ACTTGA_L002_R1_complete_filtered.fastq.gz.bigWig" )

#H3K4me3
bw_4me3 <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150366_DMSO_H3K4me3_ATCACG_L001_R1_complete_filtered.fastq.gz.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150382_SB939_H3K4me3_ACTTGA_L001_R1_complete_filtered.fastq.gz.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150350_DAC_H3K4me3_GCCAAT_L001_R1_complete_filtered.fastq.gz.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150358_DACSB_H3K4me3_CTTGTA_L001_R1_complete_filtered.fastq.gz.bigWig")

#H3K9ac
bw_9ac <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150367_DMSO_H3K9ac_CGATGT_L005_R1_complete_filtered.fastq.gz.bigWig","/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150383_SB939_H3K9ac_TGACCA_L005_R1_complete_filtered.fastq.gz.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150351_DAC_H3K9ac_TTAGGC_L005_R1_complete_filtered.fastq.gz.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150359_DACSB_H3K9ac_GATCAG_L005_R1_complete_filtered.fastq.gz.bigWig")
 
#H3K27ac
bw_27ac <- list("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150362_DMSO_H3K27ac_ATCACG_L008_R1_complete_filtered.fastq.gz.bigWig", "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150378_SB939_H3K27ac_ACTTGA_L008_R1_complete_filtered.fastq.gz.bigWig",
    "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150346_DAC_H3K27ac_GCCAAT_L008_R1_complete_filtered.fastq.gz.bigWig","/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150354_DACSB_H3K27ac_CTTGTA_L008_R1_complete_filtered.fastq.gz.bigWig")

#H3K14ac
chip.dir <- "/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/"
bw_14ac <- list(file.path(chip.dir, paste0("GSM2597665_DMSO-7.normalized", ".bigWig")),file.path(chip.dir, paste0("GSM2597647_SB-7.normalized", ".bigWig")),
    file.path(chip.dir, paste0("GSM2597638_DAC-7.normalized", ".bigWig")),file.path(chip.dir, paste0("GSM2597656_DAC+SB-7.normalized", ".bigWig")) )

#H2BK5ac
bw_5ac <- list(file.path(chip.dir, paste0("GSM2597670_DMSO-15.normalized", ".bigWig")),file.path(chip.dir, paste0("GSM2597652_SB-15.normalized", ".bigWig")),
    file.path(chip.dir, paste0("GSM2597643_DAC-15.normalized", ".bigWig")),file.path(chip.dir, paste0("GSM2597661_DAC+SB-15.normalized", ".bigWig")))

#combine lists
bw_list <- list(bw_rna, bw_cage, bw_9me3, bw_4me3, bw_9ac, bw_27ac, bw_14ac, bw_5ac)
names(bw_list) <- c("RNAseq","CAGE", "H3K9me3", "H3K4me3","H3K9ac", "H3K27ac", "H3K14ac", "H2BK5ac")

#prepare lists and sample sheets
sample_anno <- data.frame(sample_name = factor(c("DMSO", "SB", "DAC", "DACSB"),levels=c("DMSO", "SB", "DAC", "DACSB")), treatment=factor(c("DMSO", "SB", "DAC", "DACSB"),levels=c("DMSO", "SB", "DAC", "DACSB")), bigwig=NA)
#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
#select region

#get regions of interest 
anno_classi$transcript_id <- anno_classi$qry_id
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
DEG_results_list_anno_sub <- lapply(DEG_results_list, function(x){
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    x <- split(x, x$class_code_simple)
    x
})
DEG_results_list_anno_sub_sub <- DEG_results_list_anno_sub$DACandSB939_vs_DMSO

#loop over different sets of interest and plot them
dir.create(file.path(PostDE.dir, "profilePlots", "data"), recursive=TRUE)
dir.create(file.path(PostDE.dir, "profilePlots", "original_GEO"), recursive=TRUE)
dir.create(file.path(PostDE.dir, "profilePlots","original_GEO_weighted"))
dir.create(file.path(PostDE.dir, "profilePlots","original_GEO_weighted2"))

anno_transcript <- anno[anno$type=="transcript",]

for(ext in c(5000,500,1500)){
    for(i in names(DEG_results_list_anno_sub_sub)){
    trans <- DEG_results_list_anno_sub_sub[[i]]$transcript_id
    tss <- resize(anno_transcript[anno_transcript$transcript_id %in% trans, ],1)
    mcols(tss) <- NULL
    export.bed(tss, file.path(PostDE.dir, "profilePlots", "data", "tss_oi.bed"))

    for(j in names(bw_list)){
        sample_anno$bigwig <- unlist(bw_list[[j]])
        #create column data
        col_data <- read_coldata(coldata=sample_anno, files_idx=3, sample_idx=1)
        #create matrix list
        regions_matrix <- extract_matrix(coldata=col_data, 
            bed=  file.path(PostDE.dir, "profilePlots", "data", "tss_oi.bed"), 
            genome=hg19, up=ext, down=ext, binSize=10, 
            nthreads=4)
    #draw profile plot 
       # pdf(file.path(PostDE.dir, "profilePlots","original_GEO_weighted", paste(i, "_", j,"_profilePlot_ext", ext,".pdf")), height=3,5, width=3.5)
       #     print(plot_profile2(mat_list=regions_matrix, condition="treatment", weightedMean=TRUE,
       #     collapse_reps=TRUE,summarizeBy="mean", ci=FALSE, color=treat_col))
       # dev.off()
        pdf(file.path(PostDE.dir, "profilePlots","original_GEO_weighted2", paste(i, "_", j,"_profilePlot_ext", ext,".pdf")), height=3,5, width=3.5)
            print(plot_profile3(mat_list=regions_matrix, condition="treatment", weightedMean=TRUE,
            collapse_reps=TRUE,summarizeBy="mean", ci=FALSE, color=treat_col))
        dev.off()
        print(paste0(j))
    }
    print(paste0(i))
}}




#plot old tinat coordinates
dir.create(file.path(PostDE.dir, "profilePlots", "originalTinats"), recursive=TRUE)
for(ext in c(5000,500,1500)){
    for(i in names(responsive_TSS)){
    tss <- resize(responsive_TSS[[i]],1)
    mcols(tss) <- NULL
    export.bed(tss, file.path(PostDE.dir, "profilePlots", "data", "tss_oi.bed"))

    for(j in names(bw_list)){
        sample_anno$bigwig <- unlist(bw_list[[j]])
        #create column data
        col_data <- read_coldata(coldata=sample_anno, files_idx=3, sample_idx=1)
        print(sample_anno)
        #create matrix list
        regions_matrix <- extract_matrix(coldata=col_data, 
            bed=  file.path(PostDE.dir, "profilePlots", "data", "tss_oi.bed"), 
            genome=hg19, up=ext, down=ext, binSize=10, 
            nthreads=4)
    #draw profile plot 
        pdf(file.path(PostDE.dir, "profilePlots","originalTinats", paste(i, "_", j,"_profilePlot_ext", ext,".pdf")), height=3,5, width=3.5)
            print(plot_profile(mat_list=regions_matrix, condition="treatment", 
            collapse_reps=TRUE,summarizeBy="mean", ci=TRUE, color=treat_col))
        dev.off()
        print(paste0(j))
    }
    print(paste0(i))
}}



# #heaatmap
# #extract matrix for plotting
# ext <- 500
# mat_list <- lapply(bw_list, function(x){
#     normalizeToMatrix(x,repeats_LTR5all, value_column = "score" , 
#     extend = ext, mean_mode = "w0", w = 10, keep = c(0, 0.99))
# })
# #common color annotation
# common_min = min(unlist(mat_list))
# common_max = max(unlist(mat_list))
# col_fun = circlize::colorRamp2(c(common_min, common_max), c("white", "red"))
# #top legend
# pal <- c("darkgray", "orange", "tomato3")
# names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")

# lgd = Legend(at = c("LTR5_Hs", "LTR5B", "LTR5A"), title = "Repeat names", 
#     type = "lines", legend_gp = gpar(col = 2:4))


# #enriched heatmap creation
# ht_list <-  Heatmap( as.vector(repeats_LTR5all$repName), col = structure(2:4, names =c("LTR5_Hs", "LTR5B", "LTR5A")), name = "Repeat names",
#     show_row_names = FALSE, width = unit(3, "mm")) +
    
#     EnrichedHeatmap(mat_list[[1]], name=names(mat_list)[[1]], column_title = "Fuentes et al. Sp gRNAs\noriginal data",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +
    
#     EnrichedHeatmap(mat_list[[2]], name=names(mat_list)[[2]], column_title = "Non-targeting Sp gRNAs\noriginal data",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[3]], name=names(mat_list)[[3]], column_title = "Fuentes et al. Sa gRNAs\noriginal data",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[4]], name=names(mat_list)[[4]], column_title = "Fuentes et al. gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[5]], name=names(mat_list)[[5]], column_title = "RepGuide loose gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means), 
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))), 
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[6]], name=names(mat_list)[[6]], column_title = "RepGuide strict gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),  
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))),
#     col =  col_fun) +

#     EnrichedHeatmap(mat_list[[7]], name=names(mat_list)[[7]], column_title = "Non-targeting gRNAs",
#     clustering_distance_rows = dist_by_closeness,
#     #top_annotation =  HeatmapAnnotation(col_means = anno_col_means),
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), ylim=c(0,10))),  
#     col =  col_fun) 

# #plot heatmap
# dir.create(file.path(analysis.dir, "bw_heatmaps", "merged_bw"),recursive=TRUE)
# pdf(file.path(analysis.dir, "bw_heatmaps","merged_bw",paste0("enrichedHeat_bwMergedLTR5_repNameStrat_noCluster_bin10_ext", ext,"_adjFuentes.pdf")), width=20, height=7)
# draw(ht_list, split =  as.vector(repeats_LTR5all$repName),cluster_rows = F,
#     ht_gap = unit(c(2, 6, 6,6, 6, 6,6), "mm"))
# dev.off() 
