library(Biostrings)
library(ggpubr)

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

#load peptides of interest
pepteide_sub_induced <- readAAStringSet( "/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACSBvsDMSO_forJens_novel.fa")

#load transcript classification
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.mergedTranscripts.gtf.tmap"))
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id<- anno_classi$qry_id

#make dataframe of orfs, transcripts....
orf_list <- data.frame(ORF=names(pepteide_sub_induced),
    transcript_id=sapply(strsplit(sapply(strsplit(names(pepteide_sub_induced), "_",fixed=TRUE), 
        function(x)paste0(x[2],"_",x[3])) , ".p",fixed=TRUE), "[",1))
orf_list <- dplyr::left_join(orf_list, anno_classi, by="transcript_id")

#split by transcript_id
orf_list_splitClass<- split(orf_list, orf_list$class_code_simple)
#summarize class code
orf_list_splitClass_stat <- lapply(orf_list_splitClass, function(x)as.data.frame(table(x$transcript_id)))
for(i in names(orf_list_splitClass_stat)){
    orf_list_splitClass_stat[[i]]$class_code_simple <-i

}
orf_list_splitClass_stat_unlist <- do.call("rbind",orf_list_splitClass_stat)
orf_list_splitClass_stat_unlist$class_code_simple<- factor(orf_list_splitClass_stat_unlist$class_code_simple, levels=c( "chimeric (novel)", "non-chimeric (novel)", "known"))

lapply(orf_list_splitClass_stat, function(x)mean(x$Freq))
lapply(orf_list_splitClass_stat, function(x)sum(x$Freq))

#plot 
pdf(file.path("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/transdecoder_default_topStrand_8aa/", "ORFs_per_transcripts_inducedDACSB.pdf"),
    height=4, width=8)
gghistogram(orf_list_splitClass_stat_unlist, "Freq", 
    fill="class_code_simple" ,palette=class_col,add="mean", xlab="Number of ORFs per transcript",ylab="Frequency",
    combine=FALSE, bins=50, facet.by="class_code_simple", rug=TRUE)+rremove("legend")
dev.off()


# ##plot as histograms
# histo_known <- gghistogram(orf_list_splitClass_stat[["known"]], "Freq", fill=class_col["known"], xlim=c(0,710), bins=10)
# histo_chimeric <- gghistogram(orf_list_splitClass_stat[["chimeric (novel)"]], "Freq", fill=class_col["chimeric (novel)"], xlim=c(0,710), bins=2)
# histo_non_chimeric <- gghistogram(orf_list_splitClass_stat[["non-chimeric (novel)"]], "Freq", fill=class_col["non-chimeric (novel)"], xlim=c(0,710), bins=2)

# pdf(file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "ORFs_per_transcripts_inducedDACSB.pdf"),
#     height=5, width=10)
# ggarrange(histo_known, histo_chimeric,histo_non_chimeric,
#                  common.legend = TRUE,ncol = 3, nrow = 1)
# dev.off()

#calculate medium orf length
orf_width_mean <- list()
orf_width_median <- list()
for(i in names(orf_list_splitClass)){
    orf_width_mean[[i]] <- mean(width(ORFs[names(ORFs) %in% orf_list_splitClass[[i]]$ORF,]))
    orf_width_median[[i]] <- median(width(ORFs[names(ORFs) %in% orf_list_splitClass[[i]]$ORF,]))

}
ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")
