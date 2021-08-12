
library(UpSetR)
library(Biostrings)
library(ggpubr)
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#create DEG subsets for analysis
DEG_results_list <- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#subsert transcripts of interest
all <- DEG_results_list$DACandSB939_vs_DMSO
LTR12 <- all[all$dist_nearest_LTR12repeat==0,]
sign_up <- all[which(all$padj < 0.01 & all$log2FoldChange>2),]
sign_down <- all[which(all$padj < 0.01 & all$log2FoldChange<(-2)),]

DEG_list <- list(all_transcripts=all,
    all_LTR12_transcripts=all_LTR12,
    sign_downregulated_transcripts=sign_down,
    sign_upregulated_transcripts=sign_up
    )
DEG_list_transcripts <- lapply(DEG_list, rownames)


#plot upset plot
dir.create(file.path(PostDE.dir,"comparison"))
dir.create(file.path(PostDE.dir, "comparison", "transcripts"))
pdf(file.path(PostDE.dir,"comparison", "transcripts", "comparison_DEG_subsets.pdf"), width=8, height=4)
upset(fromList(DEG_list_transcripts), order.by = "freq")
dev.off()

#subset peptides based on this list
#get Peptide
peptide <- Biostrings::readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs.pep")
#remove orfs which are 3 prime or internal
length(peptide)
peptide_sub <- peptide[grep("internal", names(peptide),invert=TRUE),]
length(peptide_sub)
peptide_sub <- peptide_sub[grep("3prime", names(peptide_sub),invert=TRUE),]
#remove squence of 5 prime prior to M
peptide_5prime <- peptide_sub[grep("5prime", names(peptide_sub),invert=FALSE),]
peptide_sub <- peptide_sub[grep("5prime", names(peptide_sub),invert=TRUE),]
peptide_5prime_M <- peptide_5prime[grep("M",peptide_5prime ),]
#peptide_5prime_M <- peptide_5prime_M[1:100,]
hits = vmatchPattern("M", peptide_5prime_M, max.mismatch=0)
hits.frame=as.data.frame(hits)
peptide_5prime_M_sub <- parallel::mcmapply(subseq, peptide_5prime_M[hits.frame$group], hits.frame$start,width(peptide_5prime_M[hits.frame$group]),mc.cores=10)
peptide_5prime_M_sub <- AAStringSet(peptide_5prime_M_sub)
length(peptide_5prime_M_sub)
peptide_sub <- c(peptide_sub, peptide_5prime_M_sub)
#remove sequences after *
peptide_sub <- AAStringSet(gsub("[*]*", "", peptide_sub))
#remove sequences smaller than 8 aa
peptide_sub <- peptide_sub[width(peptide_sub)>7,]

#look at overlap
temp <- sapply(strsplit(names(peptide_sub), ".p",fixed=TRUE), "[",1) 
DEG_list_peptides <-lapply(DEG_list_transcripts, function(x){
     x <- peptide_sub[temp %in% x,]
     x 
}
)
DEG_list_peptides_names <-lapply(DEG_list_peptides,names)
lapply(DEG_list_peptides_names,length)

dir.create(file.path(PostDE.dir,"comparison"))
dir.create(file.path(PostDE.dir, "comparison", "ORFs"))
pdf(file.path(PostDE.dir,"comparison", "ORFs", "comparison_ORFs_subsets.pdf"), width=8, height=4)
upset(fromList(DEG_list_peptides_names), order.by = "freq")
dev.off()


#now get  complete union tracks for histograms
all <- DEG_results_list$DACandSB939_vs_DMSO
all_LTR12 <- all[all$dist_nearest_LTR12repeat==0,]
sign_up <- all[which(all$padj < 0.01 & all$log2FoldChange>2),]
sign_down <- all[which(all$padj < 0.01 & all$log2FoldChange<(-2)),]
sign_up_LTR12 <- sign_up[sign_up$dist_nearest_LTR12repeat==0,]
 
DEG_list <- list(all_transcripts=all,
    all_LTR12_transcripts=all_LTR12,
    sign_downregulated_transcripts=sign_down,
    sign_upregulated_transcripts=sign_up,
    sign_upregulated_LTR12_transcripts=sign_up_LTR12
    ) 

#look at orfs
temp <- sapply(strsplit(names(peptide_sub), ".p",fixed=TRUE), "[",1) 
DEG_list_peptides <-lapply(DEG_list, function(x){
     x <- peptide_sub[temp %in% rownames(x),]
     x 
}
) 

#look at length distribution
for(i in names(DEG_list_peptides)){
    size <- width(DEG_list_peptides[[i]])
    pdf(file.path(PostDE.dir,"comparison", "ORFs", paste0("histo_aminoacidLength_", i, ".pdf")))
    print(gghistogram(size, bins=50, xlab="ORF length [aminoacids]", add="mean", fill="gray", add_density = TRUE,
    ylab="Number of ORFs" ))
    dev.off()
} 
for(i in names(DEG_list_peptides)){
    size <- width(DEG_list_peptides[[i]])
    pdf(file.path(PostDE.dir,"comparison", "ORFs", paste0("histo_aminoacidLength_", i, "_log2Scale.pdf")))
    print(gghistogram(size, bins=50, xlab="ORF length [aminoacids]", add="mean", fill="gray", add_density = F,
    ylab="Number of ORFs" )+ yscale("log2", .format = TRUE))
    dev.off()
} 

#plot width vs number of orfs
size <-list()
norm_factor <- list()
for(i in names(DEG_list_peptides)){
    size[[i]] <- width(DEG_list_peptides[[i]])
    size[[i]] <- as.data.frame(table(size[[i]]))
    size[[i]]$Var1 <- as.numeric(size[[i]]$Var1)
    size[[i]]$set <- i
    pdf(file.path(PostDE.dir,"comparison", "ORFs", paste0("scatter_aminoacidLength_", i, "_log2Scale.pdf")))
    print(ggscatter(size[[i]], x="Var1", y="Freq", ylab="Number of ORFs",xlab="ORF length [aminoacids]")+ xscale("log2", .format = TRUE))
    dev.off()
}
size_ul <-do.call("rbind", size)
pdf(file.path(PostDE.dir,"comparison", "ORFs", paste0("scatter_aminoacidLength_", "combined", ".pdf")),width=10)
ggscatter(size_ul, x="Var1", y="Freq", col="set",ylab="Number of ORFs",xlab="ORF length [aminoacids]", legend="right")
dev.off()
pdf(file.path(PostDE.dir,"comparison", "ORFs", paste0("scatter_aminoacidLength_", "combined", "_log2Scale.pdf")),width=10)
ggscatter(size_ul, x="Var1", y="Freq", col="set",ylab="Number of ORFs",xlab="ORF length [aminoacids]", legend="right")+ xscale("log2", .format = TRUE)
dev.off()
#normalize size
norm_factor <- list()
for(i in names(DEG_list_peptides)){
    size[[i]] <- width(DEG_list_peptides[[i]])
    size[[i]] <- as.data.frame(table(size[[i]]))
    size[[i]]$Var1 <- as.numeric(size[[i]]$Var1)
    size[[i]]$set <- i
    norm_factor[[i]] <- sum(size[[i]]$Freq)
    size[[i]]$Freq_norm <-  size[[i]]$Freq/ norm_factor[[i]]
}
size_ul <-do.call("rbind", size)
pdf(file.path(PostDE.dir,"comparison", "ORFs", paste0("scatter_aminoacidLength_", "combined_normalized", ".pdf")),width=10)
ggscatter(size_ul, x="Var1", y="Freq_norm", col="set",ylab="Normalized number of ORFs",xlab="ORF length [aminoacids]", legend="right")
dev.off()

pdf(file.path(PostDE.dir,"comparison", "ORFs", paste0("scatter_aminoacidLength_", "combined_normalized", "_log2Scale.pdf")),width=10)
ggscatter(size_ul, x="Var1", y="Freq_norm", col="set",ylab="Normalized number of ORFs",xlab="ORF length [aminoacids]", legend="right")+ xscale("log2", .format = TRUE)
dev.off()

