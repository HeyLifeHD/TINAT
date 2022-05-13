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
pepteide_sub_induced <- readAAStringSet( "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")

#load transcript classification
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
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
mean(do.call("rbind", orf_list_splitClass_stat)$Freq)
lapply(orf_list_splitClass_stat, function(x)sum(x$Freq))

#plot 
pdf(file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "ORFs_per_transcripts_inducedDACSB.pdf"),
    height=4, width=8)
gghistogram(orf_list_splitClass_stat_unlist, "Freq", 
    fill="class_code_simple" ,palette=class_col,add="mean", xlab="Number of ORFs per transcript",ylab="Frequency",
    combine=FALSE, bins=50, facet.by="class_code_simple", rug=TRUE)+rremove("legend")
dev.off()

pdf(file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "ORFs_per_transcripts_inducedDACSB_comb.pdf"),
    height=4, width=4)
gghistogram(orf_list_splitClass_stat_unlist, "Freq", 
    color="class_code_simple" ,palette=class_col,add="mean", xlab="Number of ORFs per transcript",ylab="Frequency",
    bins=50, rug=TRUE)+rremove("legend")
dev.off()

pdf(file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "ORFs_per_transcripts_inducedDACSB_comb_noCol.pdf"),
    height=4, width=4)
gghistogram(orf_list_splitClass_stat_unlist, "Freq", 
    add="mean", xlab="Number of ORFs per transcript",ylab="Frequency",
    combine=FALSE, bins=50, rug=TRUE)+rremove("legend")
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



#plot length vs number of orfs
ORFs <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")
ORFs_uniprot <- readAAStringSet("/omics/groups/OE0219/internal/tinat/integration/peptidomics/data/190821_uniprot_homo sapiens.fasta")
#make dataframe of orfs, transcripts....
orf_list <- data.frame(ORF=names(ORFs),ORFLength=width(ORFs),
    transcript_id=sapply(strsplit(sapply(strsplit(names(ORFs), "_",fixed=TRUE), 
        function(x)paste0(x[2],"_",x[3])) , ".p",fixed=TRUE), "[",1), class="DIT ORFs")
orf_list <- dplyr::left_join(orf_list, anno_classi, by="transcript_id")
orf_uniprot_list <- data.frame(ORF=names(ORFs_uniprot),ORFLength=width(ORFs_uniprot),class="Uniprot ORFs")
#get stats
stat <- as.data.frame(table(orf_list$ORFLength))
stat$Var1 <- as.integer(levels(stat$Var1))
stat$class <- "DIT ORFs"
stat_uniprot <- as.data.frame(table(orf_uniprot_list$ORFLength))
stat_uniprot$Var1 <- as.integer(levels(stat_uniprot$Var1))
stat_uniprot$class <- "Uniprot ORFs"
stat_combined <- rbind(stat, stat_uniprot)
fig <- ggscatter(stat_combined, x="Var1", y="Freq", xlab="ORF length [AA]",col="class",
    fill="class",size=0.3, palette=c("gray", "black"),ylab=  "Number of ORFs")+
    rremove("legend")+
    #coord_trans(x="log10")
    scale_x_continuous(trans='log10')
ggsave(plot=fig, filename=file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "Number_vs_length_ORFs_all.pdf"),useDingbats=FALSE, height=4, width=4)
fig <- ggscatter(stat_combined, x="Var1", y="Freq", xlab="ORF length [AA]",col="class",
    fill="class",size=0.3, palette=c("gray", "black"),ylab=  "Number of ORFs")+
    rremove("legend")+
    #coord_trans(x="log10")
    scale_x_continuous(trans='log10')+
    scale_y_continuous(trans='log10')

ggsave(plot=fig, filename=file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "Number_vs_length_ORFs_all_V2.pdf"),useDingbats=FALSE, height=4, width=4)

#stratified by class code
orf_list_split <- split(orf_list, orf_list$class_code_simple)
stat <- lapply(names(orf_list_split), function(x){
    stat <- as.data.frame(table(orf_list_split[[x]]$ORFLength))
    stat$Var1 <- as.integer(levels(stat$Var1))
    stat$class_code <- x
    stat
})
stat <- do.call("rbind",stat)
fig <- ggscatter(stat, x="Var1", y="Freq", xlab="ORF length [AA]",col="class_code",facet.by="class_code",fill="class_code",size=0.3, palette=class_col,ylab=  "Number of ORFs")+rremove("legend")+scale_x_continuous(trans='log10')
ggsave(plot=fig, filename=file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "Number_vs_length_ORFs_stratClassCode.pdf"),useDingbats=FALSE, height=4, width=12)
fig <- ggscatter(stat, x="Var1", y="Freq", xlab="ORF length [AA]",col="class_code",facet.by="class_code",fill="class_code",size=0.3, palette=class_col,ylab=  "Number of ORFs")+rremove("legend")+scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')
ggsave(plot=fig, filename=file.path("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/", "Number_vs_length_ORFs_stratClassCode_v2.pdf"),useDingbats=FALSE, height=4, width=12)
