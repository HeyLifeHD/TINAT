
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

#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 4##set logfold2 cutoff


#subsert transcripts of interest
all <- DEG_results_list$DACandSB939_vs_DMSO
LTR12 <- all[all$dist_nearest_LTR12repeat==0,]
sign_up <- all[which(all$padj < alpha & all$log2FoldChange>lfc),]
sign_down <- all[which(all$padj < alpha & all$log2FoldChange<(-lfc)),]

DEG_list <- list(all_transcripts=all,
    all_LTR12_transcripts=all_LTR12,
    sign_downregulated_transcripts=sign_down,
    sign_upregulated_transcripts=sign_up
    )
DEG_list_transcripts <- lapply(DEG_list, rownames)


#plot upset plot
dir.create(file.path(PostDE.dir,"comparison"))
dir.create(file.path(PostDE.dir,"comparison","4lfc"))
dir.create(file.path(PostDE.dir, "comparison","4lfc", "transcripts"))
pdf(file.path(PostDE.dir,"comparison", "4lfc", "transcripts", "comparison_DEG_subsets.pdf"), width=8, height=4)
upset(fromList(DEG_list_transcripts), order.by = "freq", set_size.show = TRUE)
dev.off()
 


#now get  complete union tracks for histograms
all <- DEG_results_list$DACandSB939_vs_DMSO
all_LTR12 <- all[all$dist_nearest_LTR12repeat==0,]
sign_up <- all[which(all$padj < alpha & all$log2FoldChange>lfc),]
sign_down <- all[which(all$padj < alpha & all$log2FoldChange<(-lfc)),]
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
    pdf(file.path(PostDE.dir,"comparison","4lfc", "ORFs", paste0("histo_aminoacidLength_", i, ".pdf")))
    print(gghistogram(size, bins=50, xlab="ORF length [aminoacids]", add="mean", fill="gray", add_density = TRUE,
    ylab="Number of ORFs" ))
    dev.off()
} 
for(i in names(DEG_list_peptides)){
    size <- width(DEG_list_peptides[[i]])
    pdf(file.path(PostDE.dir,"comparison","4lfc", "ORFs", paste0("histo_aminoacidLength_", i, "_log2Scale.pdf")))
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
    pdf(file.path(PostDE.dir,"comparison","4lfc", "ORFs", paste0("scatter_aminoacidLength_", i, "_log2Scale.pdf")))
    print(ggscatter(size[[i]], x="Var1", y="Freq", ylab="Number of ORFs",xlab="ORF length [aminoacids]")+ xscale("log2", .format = TRUE))
    dev.off()
}
size_ul <-do.call("rbind", size)
pdf(file.path(PostDE.dir,"comparison",,"4lfc" "ORFs", paste0("scatter_aminoacidLength_", "combined", ".pdf")),width=10)
ggscatter(size_ul, x="Var1", y="Freq", col="set",ylab="Number of ORFs",xlab="ORF length [aminoacids]", legend="right")
dev.off()
pdf(file.path(PostDE.dir,"comparison","4lfc", "ORFs", paste0("scatter_aminoacidLength_", "combined", "_log2Scale.pdf")),width=10)
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
pdf(file.path(PostDE.dir,"comparison","4lfc", "ORFs", paste0("scatter_aminoacidLength_", "combined_normalized", ".pdf")),width=10)
ggscatter(size_ul, x="Var1", y="Freq_norm", col="set",ylab="Normalized number of ORFs",xlab="ORF length [aminoacids]", legend="right")
dev.off()

pdf(file.path(PostDE.dir,"comparison","4lfc", "ORFs", paste0("scatter_aminoacidLength_", "combined_normalized", "_log2Scale.pdf")),width=10)
ggscatter(size_ul, x="Var1", y="Freq_norm", col="set",ylab="Normalized number of ORFs",xlab="ORF length [aminoacids]", legend="right")+ xscale("log2", .format = TRUE)
dev.off()

