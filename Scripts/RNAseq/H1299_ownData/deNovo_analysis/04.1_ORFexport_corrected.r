#change orf naming for Jens gusto
library(Biostrings)
peptide <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs.pep")

#get unique 
unique(lapply(strsplit(names(peptide), " ",fixed=TRUE), "[", 2))

#get complete orfs
peptide_sub <- peptide[grep("complete", names(peptide),invert=FALSE),]
length(peptide_sub)

#get 5 prime orfs
peptide_5prime <- peptide[grep("5prime", names(peptide),invert=FALSE),]
length(peptide_5prime)
#find M
peptide_5prime_M <- peptide_5prime[grep("M",peptide_5prime ),]
length(peptide_5prime_M)
hits <-  vmatchPattern("M", peptide_5prime_M, max.mismatch=0)
#only keep first hit
hits_first <- mclapply(hits, function(x)x[1,],mc.cores=10)
length(hits_first)==length(peptide_5prime_M)
#remove squence of 5 prime prior to M
hits_frame <- unlist(IRangesList(hits_first))
length(hits_frame)==length(peptide_5prime_M)
hits_frame <- as.data.frame(hits_frame)
peptide_5prime_M_sub <- parallel::mcmapply(subseq, peptide_5prime_M[hits_frame$names], hits_frame$start,width(peptide_5prime_M[hits_frame$names]),mc.cores=10)
peptide_5prime_M_sub <- AAStringSet(peptide_5prime_M_sub)
length(peptide_5prime_M_sub)==length(peptide_5prime_M)

#combine sequences 
peptide_sub <- c(peptide_sub, peptide_5prime_M_sub)

#remove sequences after *
peptide_sub <- AAStringSet(gsub("[*]*", "", peptide_sub))

#remove sequences smaller than 8 aa
peptide_sub <- peptide_sub[width(peptide_sub)>7,]
subset(peptide_sub, 1)
#QC
length(peptide_sub)
table(duplicated(names(peptide_sub)))
table(subseq(peptide_sub, start = 1, end = 1)=="M")

#rename orfs
names(peptide_sub) <- paste0("ORF_", names(peptide_sub),"_tinat")
names(peptide_sub)  <- gsub("ORF_ORF_", "ORF_", names(peptide_sub))
names(peptide_sub)  <- gsub("_tinat_tinat", "_tinat", names(peptide_sub))

writeXStringSet(peptide_sub, "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_forJens_novel.fa")


##select ORFs baased on upregulation in DEG analysis
DEG_results_list<- readRDS("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/results/PostDE/DEG_results_group_list.rds")
induced_genes <- DEG_results_list$DACandSB939_vs_DMSO[which(DEG_results_list$DACandSB939_vs_DMSO$padj < 0.01 & DEG_results_list$DACandSB939_vs_DMSO$log2FoldChange>2),]

#subset peptide list
temp <- sapply(strsplit(names(peptide_sub), "_",fixed=TRUE), function(x)paste0(x[2],"_",x[3])) 
temp <- sapply(strsplit(temp, ".p",fixed=TRUE), "[",1) 
pepteide_sub_induced <- peptide_sub[temp %in% rownames(induced_genes),]
length(pepteide_sub_induced)
length(peptide_sub)

writeXStringSet(pepteide_sub_induced, "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")


##select ORFs baased on novel and upregulation in DEG analysis
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)


induced_genes <- induced_genes[rownames(induced_genes) %in% anno_classi[anno_classi$class_code_simple !="known",]$qry_id, ]

#subset peptide list
#temp <- sapply(strsplit(names(peptide_sub), "_",fixed=TRUE), function(x)paste0(x[2],"_",x[3])) 
#temp <- sapply(strsplit(temp, ".p",fixed=TRUE), "[",1) 
pepteide_sub_induced <- peptide_sub[temp %in% rownames(induced_genes),]
length(pepteide_sub_induced)
length(peptide_sub)
length(peptide)

writeXStringSet(pepteide_sub_induced, "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_novelTranscript_forJens_novel.fa")


