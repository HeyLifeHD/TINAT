#libraries
library(Biostrings)

#load data
ORFs_cellline <- readAAStringSet("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACSBvsDMSO_forJens_novel.fa")
ORFs_H1299 <- readAAStringSet("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")

#rename orfs
names(ORFs_cellline)<- sapply(strsplit(names(ORFs_cellline), " type", fixed=TRUE), "[",1)
names(ORFs_H1299)<- sapply(strsplit(names(ORFs_H1299), " type", fixed=TRUE), "[",1)

#subset orfs of interest
set_cellline <-c("ORF_MSTRG.3332.1.p1", "ORF_MSTRG.3361.1.p1", "ORF_MSTRG.11293.1.p12", "ORF_MSTRG.11293.2.p12", "ORF_MSTRG.11293.3.p12", "ORF_MSTRG.11293.4.p13", "ORF_MSTRG.11293.5.p13", "ORF_MSTRG.11293.6.p13", "ORF_MSTRG.21956.1.p3", "ORF_MSTRG.21956.2.p1", "ORF_MSTRG.21956.5.p2", "ORF_MSTRG.23950.1.p19", "ORF_MSTRG.23950.2.p1", "ORF_MSTRG.23950.2.p20", "ORF_MSTRG.25549.1.p1", "ORF_MSTRG.30823.1.p1", "ORF_MSTRG.31089.1.p2", "ORF_MSTRG.32391.3.p4", "ORF_MSTRG.32391.4.p4", "ORF_MSTRG.32391.5.p5", "ORF_MSTRG.32528.1.p3", "ORF_MSTRG.32528.2.p4", "ORF_MSTRG.32553.1.p1", "ORF_MSTRG.33292.4.p2", "ORF_MSTRG.33292.1.p1", "ORF_MSTRG.33292.2.p1", "ORF_MSTRG.33292.3.p1", "ORF_MSTRG.33292.6.p1", "ORF_MSTRG.33292.8.p1", "ORF_MSTRG.33292.10.p1", "ORF_MSTRG.33292.14.p2", "ORF_ENST00000500169.6_1.p2", "ORF_MSTRG.34228.1.p1", "ORF_MSTRG.34248.2.p1", "ORF_MSTRG.34463.1.p1", "ORF_MSTRG.39579.3.p1", "ORF_MSTRG.39841.1.p2", "ORF_MSTRG.39845.1.p1", "ORF_MSTRG.39845.16.p1", "ORF_MSTRG.39906.8.p2", "ORF_MSTRG.39906.12.p2", "ORF_MSTRG.42343.1.p15", "ORF_MSTRG.42380.1.p2", "ORF_MSTRG.42380.2.p2", "ORF_MSTRG.42379.1.p1", "ORF_MSTRG.44444.1.p6", "ORF_MSTRG.45990.1.p1", "ORF_MSTRG.45990.2.p1", "ORF_MSTRG.45990.3.p1", "ORF_MSTRG.45990.4.p1", "ORF_MSTRG.45990.6.p1", "ORF_MSTRG.45990.7.p1", "ORF_MSTRG.45990.9.p1", "ORF_MSTRG.45990.10.p1", "ORF_MSTRG.46862.1.p5")
length(set_cellline)
set_h1299 <- c("ORF_MSTRG.30117.3.p1", "ORF_ENST00000500169.6_1.p1", "ORF_MSTRG.24557.3.p4", "ORF_MSTRG.2576.1.p1", "ORF_MSTRG.18454.1.p1", "ORF_MSTRG.25979.1.p1", "ORF_MSTRG.8418.1.p10", "ORF_MSTRG.16894.1.p1", "ORF_MSTRG.16894.3.p1", "ORF_MSTRG.18454.1.p2")
length(set_h1299)

#subset ORFs of interest
ORFs_cellline_sub <- ORFs_cellline[set_cellline,]
ORFs_H1299_sub <- ORFs_H1299[set_h1299,]
ORF_list <- list(H1299=ORFs_H1299_sub, Cellline=ORFs_cellline_sub)

#look for peptides
sum(vcountPattern("TPSIRTVTL",ORFs_cellline))
sum(vcountPattern("TPSLRTVTL",ORFs_cellline))
sum(vcountPattern("TPSIRTVTL",ORFs_H1299))
sum(vcountPattern("TPSLRTVTL",ORFs_H1299))
allH1299_leucine_match <- ORFs_H1299[vcountPattern("TPSLRTVTL",ORFs_H1299)>0,]
paste0(names(allH1299_leucine_match), collapse=", ")
paste0(names(allH1299_leucine_match)[!names(allH1299_leucine_match) %in% c("ORF_MSTRG.30117.3.p1", "ORF_ENST00000500169.6_1.p1", "ORF_MSTRG.24557.3.p4", "ORF_MSTRG.2576.1.p1", "ORF_MSTRG.18454.1.p1", "ORF_MSTRG.25979.1.p1", "ORF_MSTRG.8418.1.p10", "ORF_MSTRG.16894.1.p1", "ORF_MSTRG.16894.3.p1", "ORF_MSTRG.18454.1.p2")], collapse=", ") 

isoleucine_count <- lapply(ORF_list, function(x){
    y <- vcountPattern("TPSIRTVTL",x)
    data.frame(ORF=names(x), Frequency=y)
})
leucine_count <- lapply(ORF_list, function(x){
    y <- vcountPattern("TPSLRTVTL",x)
    data.frame(ORF=names(x), Frequency=y)
})
H1299_counts <- cbind(isoleucine_count$H1299,leucine_count$H1299)[,c(1,2,4)]
colnames(H1299_counts)<- c("ORF", "IsoleucineFrequency", "LeucineFrequency")
H1299_counts$ORF <- as.character(H1299_counts$ORF)
H1299_counts <- rbind(H1299_counts, c(c(ORF="Total"), colSums(H1299_counts[,c(2,3)])))
cellline_counts <- cbind(isoleucine_count$Cellline,leucine_count$Cellline)[,c(1,2,4)]
colnames(cellline_counts)<- c("ORF", "IsoleucineFrequency", "LeucineFrequency")
cellline_counts$ORF <- as.character(cellline_counts$ORF)
cellline_counts <- rbind(cellline_counts, c(c(ORF="Total"), colSums(cellline_counts[,c(2,3)])))

dir.create("/omics/groups/OE0219/internal/tinat/temp")
write.table(H1299_counts,"/omics/groups/OE0219/internal/tinat/temp/H1299_counts.tsv", sep="\t", row.names=FALSE, col.names=TRUE,quote=FALSE)
write.table(cellline_counts,"/omics/groups/OE0219/internal/tinat/temp/cellline_counts.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)