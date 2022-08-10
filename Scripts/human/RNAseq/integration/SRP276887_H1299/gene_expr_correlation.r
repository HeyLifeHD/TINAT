#comparison of H1299 with setb1 data
#libraries
library(ggpubr)

#directories
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/H1299_SRP276887/"
dir.create(output.dir, recursive=TRUE)

#load results
DEG_results_list_H1299 <- readRDS("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/results/PostDE/DEG_results_group_list.rds")
DEG_results_list_SETB1 <- readRDS("/omics/groups/OE0219/internal/tinat/220118_SRP276887_deNovo_quantification_H12Ref_analysis/results/PostDE/DEG_results_group_list.rds")
#load anno classification
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))


#combine data frames of main comparison
merged_results <- dplyr::left_join(DEG_results_list_H1299$DACandSB939_vs_DMSO,DEG_results_list_SETB1$SETB1_KO_vs_Control, by="transcript_id" )
#merge with annotation
anno_classi_sub <- anno_classi[, c("class_code","qry_id", "num_exons" )]
colnames(anno_classi_sub)<- c("class_code", "transcript_id", "num_exons")
merged_results <- dplyr::left_join(merged_results, anno_classi_sub)

#plot correlation
pal <- randomcoloR::distinctColorPalette(length(unique(merged_results$class_code)))
names(pal)<- unique(merged_results$class_code)
#pearson

g1 <- ggscatter(merged_results, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",size=.3,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control",
            add = "reg.line",conf.int = T)  + guides(color="none") +  
        stat_cor(label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_all_genes.pdf"), height = 6, width = 7)
print(g1)
dev.off()

#plot correlation colored by class anno
g1 <- ggscatter(merged_results, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",color = "class_code",palette=pal,size=.3,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control",
            add = "reg.line",conf.int = T)  + guides(color="none") +  
        stat_cor(aes(color=class_code),label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_all_genes_classCode.pdf"), height = 6, width = 7)
print(g1)
dev.off()

#plot correlation colored by class anno only novel
merged_results_sub <- merged_results[merged_results$class_code!= "=",]
g1 <- ggscatter(merged_results_sub, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",color = "class_code",palette=pal,size=.3,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control",
            add = "reg.line",conf.int = T)  + guides(color="none") +  
        stat_cor(aes(color=class_code),label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_all_genes_classCode_novel.pdf"), height = 6, width = 7)
print(g1)
dev.off()

#subset genes significantly upregulated in DAC_SB
merged_results_sub <- merged_results[which(merged_results$log2FoldChange.x >2 & merged_results$padj.x < 0.01),]
g1 <- ggscatter(merged_results_sub, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",color = "class_code",palette=pal,size=.3,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control",
            add = "reg.line",conf.int = T)  + guides(color="none") +  
        stat_cor(aes(color=class_code),label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_DACSB_DEGup_classCode.pdf"), height = 6, width = 7)
print(g1)
dev.off()
#and novel transcript
merged_results_sub <- merged_results[which(merged_results$log2FoldChange.x >2 & merged_results$padj.x < 0.01),]
merged_results_sub <- merged_results_sub[merged_results_sub$class_code!= "=",]

g1 <- ggscatter(merged_results_sub, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",color = "class_code",palette=pal,size=.3,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control",
            add = "reg.line",conf.int = T)  + guides(color="none") +  
        stat_cor(aes(color=class_code),label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_DACSB_DEGup_classCode_novel.pdf"), height = 6, width = 7)
print(g1)
dev.off()

# novel and LTR derived transcript
merged_results_sub <- merged_results[which(merged_results$log2FoldChange.x >2 & merged_results$padj.x < 0.01),]
merged_results_sub <- merged_results_sub[merged_results_sub$class_code!= "=",]
merged_results_sub <- merged_results_sub[merged_results_sub$dist_nearest_LTR12repeat.x == 0,]
dim(merged_results_sub)

g1 <- ggscatter(merged_results_sub, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",color = "class_code",palette=pal,size=.3,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control",
            add = "reg.line",conf.int = T)  + guides(color="none") +  
        stat_cor(aes(color=class_code),label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_DACSB_DEGup_classCode_novel_LTR12derived.pdf"), height = 6, width = 7)
print(g1)
dev.off()


#subset genes significantly upregulated in Setb1
merged_results_sub <- merged_results[which(merged_results$log2FoldChange.y >1 & merged_results$padj.y < 0.05),]
merged_results_sub <- merged_results_sub[!is.na(merged_results_sub$log2FoldChange.x)|is.na(merged_results_sub$log2FoldChange.y),]
g1 <- ggscatter(merged_results_sub, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",color = "class_code",palette=pal,size=.3,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control",
            add = "reg.line",conf.int = T)  + guides(color="none") +  
        stat_cor(aes(color=class_code),label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_SETB1_DEGup_classCode.pdf"), height = 6, width = 7)
print(g1)
dev.off()



#subset genes significantly upregulated in Setb1 and LTR12 derived
merged_results_sub <- merged_results[which(merged_results$log2FoldChange.y >1 & merged_results$padj.y < 0.05),]
merged_results_sub <- merged_results_sub[!is.na(merged_results_sub$log2FoldChange.x)|is.na(merged_results_sub$log2FoldChange.y),]
merged_results_sub <- merged_results_sub[merged_results_sub$dist_nearest_LTR12repeat.x == 0,]
dim(merged_results_sub)
g1 <- ggscatter(merged_results_sub, x="log2FoldChange.x", y="log2FoldChange.y",legend="right",color = "class_code",palette=pal,size=1,
            xlab="Gene expression change [log2FC]\n DAC+SB939 vs. DMSO" , ylab="Gene expression change [log2FC]\n SETB1-KO vs. control", #add = "reg.line"
           ,conf.int = F)  + guides(color="none") +  
        stat_cor(aes(color=class_code),label.y.npc = "bottom", label.x.npc = .6) # + yscale("log2", .format = TRUE)+xscale("log2", .format = TRUE)
pdf(file.path(output.dir, "correlation_SETB1_DEGup_classCode_LTR12derived.pdf"), height = 6, width = 7)
print(g1)
dev.off()
