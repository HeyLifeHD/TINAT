#change orf naming for Jens gusto
library(Biostrings)
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

#rename orfs
names(peptide_sub) <- paste0("ORF_", names(peptide_sub),"_tinat")

writeXStringSet(peptide_sub, "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_forJens.fa")


##select ORFs baased on upregulation in DEG analysis
DEG_results_list<- readRDS("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/results/PostDE/DEG_results_group_list.rds")
induced_genes <- DEG_results_list$DACandSB939_vs_DMSO[which(DEG_results_list$DACandSB939_vs_DMSO$padj < 0.01 & DEG_results_list$DACandSB939_vs_DMSO$log2FoldChange>2),]

#subset peptide list
temp <- sapply(strsplit(names(peptide_sub), "_",fixed=TRUE), function(x)paste0(x[2],"_",x[3])) 
temp <- sapply(strsplit(temp, ".p",fixed=TRUE), "[",1) 
pepteide_sub_induced <- peptide_sub[temp %in% rownames(induced_genes),]
length(pepteide_sub_induced)
length(peptide_sub)
length(peptide)

writeXStringSet(pepteide_sub_induced, "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens.fa")



##select ORFs baased on upregulation in DEG analysis as well as overlap to LTR12 element
DEG_results_list<- readRDS("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/results/PostDE/DEG_results_group_list.rds")
induced_genes_LTR <- induced_genes[induced_genes$dist_nearest_LTR12repeat==0,]
dim(induced_genes)
dim(induced_genes_LTR)

#subset peptide list
temp <- sapply(strsplit(names(peptide_sub), "_",fixed=TRUE), function(x)paste0(x[2],"_",x[3])) 
temp <- sapply(strsplit(temp, ".p",fixed=TRUE), "[",1) 
pepteide_sub_induced_LTR <- peptide_sub[temp %in% rownames(induced_genes_LTR),]
length(pepteide_sub_induced_LTR)
length(pepteide_sub_induced)
length(peptide_sub)
length(peptide)

writeXStringSet(pepteide_sub_induced_LTR, "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_LTR12TSS_forJens.fa")
