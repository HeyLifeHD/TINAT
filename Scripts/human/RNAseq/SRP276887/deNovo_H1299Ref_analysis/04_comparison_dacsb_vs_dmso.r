
#Directories
#setb1
base_setb1.dir<- "/omics/groups/OE0219/internal/tinat/220118_SRP276887_deNovo_quantification_H12Ref_analysis/"
base_results_setb1.dir <- file.path(base_setb1.dir, "results")
results_setb1.dir<- file.path(base_results_setb1.dir , "tables")
PreDE_setb1.dir <- file.path(base_results_setb1.dir,"PreDE")
PostDE_setb1.dir <- file.path(base_results_setb1.dir,"PostDE")
#h1299
base_h1299.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results_h1299.dir <- file.path(base_h1299.dir, "results")
results_h1299.dir<- file.path(base_results_h1299.dir , "tables")
PreDE_h1299.dir <- file.path(base_results_h1299.dir,"PreDE")
PostDE_h1299.dir <- file.path(base_results_h1299.dir,"PostDE")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
#subset anno
anno_transcript <- anno[anno$type=="transcript",]
anno_transcript <- as.data.frame(anno_transcript)
#Read in Data
DEG_results_list_setb1<- readRDS(file.path(PostDE_setb1.dir, "DEG_results_group_list.rds"))
DEG_results_list_h1299<- readRDS(file.path(PostDE_h1299.dir, "DEG_results_group_list.rds"))

#subset induced genes
induced_setb1 <- DEG_results_list_setb1$SETB1_KO_vs_Control[which(DEG_results_list_setb1$SETB1_KO_vs_Control$padj < 0.05 & 
    DEG_results_list_setb1$SETB1_KO_vs_Control$log2FoldChange>1),]
induced_setb1 <- dplyr::left_join(induced_setb1, anno_transcript, by="transcript_id")
write.table(induced_setb1, file.path(PostDE_setb1.dir,"SETB1_DEGs_induced_in_DACSB.txt"), row.names=FALSE, col.names=TRUE,sep="\t" )
induced_h1299 <- DEG_results_list_h1299$DACandSB939_vs_DMSO[which(DEG_results_list_h1299$DACandSB939_vs_DMSO$padj < 0.01 & 
    DEG_results_list_h1299$DACandSB939_vs_DMSO$log2FoldChange>2),]
induced_h1299 <- dplyr::left_join(induced_h1299, anno_transcript, by="transcript_id")
write.table(induced_h1299, file.path(PostDE_setb1.dir,"DACSB_DEGs_induced_in_SETB1.txt"), row.names=FALSE, col.names=TRUE,sep="\t" )

table(rownames(induced_setb1) %in% rownames(induced_h1299))


pdf(file.path(PostDE_setb1.dir, "Venn_DEG_up_DACSB.pdf"))
ggVennDiagram::ggVennDiagram(list(SETB1_induced=rownames(induced_setb1), DACSB_induced=rownames(induced_h1299) ), label_alpha = 0)
    #scale_fill_gradient(low="white",high = "#EF8A62")
dev.off()