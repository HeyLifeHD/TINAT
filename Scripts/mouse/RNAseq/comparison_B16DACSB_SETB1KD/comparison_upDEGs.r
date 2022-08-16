#Perform differential epression analysis of different contrasts
#libraries
library
#Directories
#own data
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_analysis"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Read in Data
DEG_results_list_DACSB<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#setb1
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220811_RNAseqSETB1KD_deNovoB16DACSB_analysis"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Read in Data
DEG_results_list_SETB1<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))




#subset upregulated genes
DEG_results_list_DACSB_up_id <- lapply(DEG_results_list_DACSB, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange>lfc),]
    x$transcript_id
})
DEG_results_list_SETB1_up_id <- lapply(DEG_results_list_SETB1, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange>lfc),]
    x$transcript_id
})
#get transcript ids of sets of interest
list_oi <- list(DACSB_vs_DMSO=DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO, SETB1KD_vs_control=DEG_results_list_SETB1_up_id$SETDB1_vs_control)

dir.create("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD", "Venn_DEG_up.pdf"))
ggVennDiagram::ggVennDiagram(list_oi, label_alpha = 0, category.names=names(list_oi))+
    rremove("legend") + 
    theme(text = element_text(size = 8)) #+ 
    #scale_fill_gradient(low="white",high = "#EF8A62")
dev.off()

write.table(DEG_results_list_DACSB$DACandSB939_vs_DMSO[DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO[DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO %in%  DEG_results_list_SETB1_up_id$SETDB1_vs_control],], 
file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD", "common_up_DACSBStats.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)