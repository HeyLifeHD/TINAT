#deg up heatmap
#libraries
library(bsseq)
library(pheatmap)
library(ggpubr)
library(readxl)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap) 

set.seed(42)

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.mergedTranscripts.gtf.tmap"))

#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACSB")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
cell_col  <- brewer.pal(length(unique(colData(vst)$cellline)), "Paired")
names(cell_col)<- as.character(unique(colData(vst)$cellline))
#anno
anno_colors <- list(treatment=treat_col, cellline= cell_col)
anno <- colData(vst)
annovst <- anno[,c("cellline", "treatment" ), drop=FALSE]


#heatmaps
#define padj cutoff for plotting --> only up
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff
cutoff <- alpha
l2fc <- lfc

#get class code information of de novo assembly, joint with number of exons
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id

DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACSB_vs_DMSO")]
genes2plot <- lapply(DEG_results_list_sub, function(x){
  x <- dplyr::left_join(x, anno_classi, by="transcript_id")
  x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
  print(dim(x))
  x <- split(x, x$class_code_simple.y)
  x
})

# complex heatmap style
#get top annotation
for(i in names(genes2plot)){
    pheno <- colData(vst)
    plot <- assay(vst)[genes2plot[[i]][[1]]$transcript_id,]
    ha <- HeatmapAnnotation(treatment = as.character(pheno$treatment),cellline= as.character(pheno$cellline), 
        col = list(treatment=treat_col[as.character(pheno$treatment)], cellline=cell_col[as.character(pheno$cellline)] ))

    #prepare heatmaps
    mat <- assay(vst)[genes2plot[[i]][["chimeric (novel)"]]$transcript_id,]
    rownames(mat) <- NULL
    ht1 <- Heatmap(mat, show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression",  cluster_columns=TRUE,show_column_names=FALSE,show_column_dend=FALSE,
            row_title="chinmeric (novel)",col = c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE,# show_row_names=FALSE, 
            meth=anno_boxplot(mat, pch=".", height= unit(1, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["chimeric (novel)"]), width=unit(.2, "cm") )))

    mat <- assay(vst)[genes2plot[[i]][["non-chimeric (novel)"]]$transcript_id,]
    rownames(mat) <- NULL
    ht2 <- Heatmap(mat, show_row_dend=FALSE,show_column_dend=FALSE, name = "z Scaled\ntranscript\nexpression",show_column_names=FALSE,
            cluster_columns=TRUE, row_title="non-chimeric (novel)",col = c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))), 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE,#show_row_names=FALSE, 
            meth=anno_boxplot(mat, pch=".", height= unit(1, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["non-chimeric (novel)"]), width=unit(.2, "cm") )))
    
     mat <- assay(vst)[genes2plot[[i]][["known"]]$transcript_id,]
    rownames(mat) <- NULL
    ht3 <- Heatmap(mat, show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression", show_column_names=FALSE,
            cluster_columns=TRUE,show_column_dend=FALSE, row_title="known",col = c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))), 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE, #show_row_names=FALSE, 
            meth=anno_boxplot(mat, pch=".", height= unit(1, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["known"]), width=unit(.2, "cm") )))

    #combine
    ht_list = ha %v% ht1 %v% ht2 %v% ht3
    #draw heatmap
    pdf(file.path(PostDE.dir, i, paste0("Complexheatmap_DEGup_stratClassCode.pdf")), height=7.5, width=5)
    draw(ht_list)
    dev.off()
    print(i)
}

#rowscaled
# complex heatmap style
#get top annotation

for(i in names(genes2plot)){
    pheno <- colData(vst)[,]
    plot <- assay(vst)[genes2plot[[i]][[1]]$transcript_id,]
    ha <- HeatmapAnnotation(treatment = as.character(pheno$treatment),cellline= as.character(pheno$cellline), 
        col = list(treatment=treat_col[as.character(pheno$treatment)], cellline=cell_col[as.character(pheno$cellline)] ))

    #prepare heatmaps
    mat <- assay(vst)[genes2plot[[i]][["chimeric (novel)"]]$transcript_id,]
    rownames(mat) <- NULL
    ht1 <- Heatmap(scale_rows(mat), show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression", cluster_columns=TRUE,show_column_names=FALSE,show_column_dend=FALSE,
            row_title="chinmeric (novel)",col =c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))) , 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE,# show_row_names=FALSE, 
            meth=anno_boxplot(mat, ylim=c(0,17), border=TRUE, outline=FALSE, pch=".", height= unit(1.5, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["chimeric (novel)"],fontsize=6,annotation_nam=gpar(fontsize=6)), width=unit(.2, "cm") )))

    mat <- assay(vst)[genes2plot[[i]][["non-chimeric (novel)"]]$transcript_id,]
    rownames(mat) <- NULL
    ht2 <- Heatmap(scale_rows(mat), show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression",cluster_columns=TRUE,show_column_names=FALSE,show_column_dend=FALSE,
            row_title="non-chimeric (novel",col = c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))) , 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE,#show_row_names=FALSE, 
            meth=anno_boxplot(mat, ylim=c(0,17), border=TRUE, outline=FALSE, pch=".", height= unit(1.5, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["non-chimeric (novel)"], fontsize=6,annotation_nam=gpar(fontsize=6)), width=unit(.2, "cm") )))
    
     mat <- assay(vst)[genes2plot[[i]][["known"]]$transcript_id,]
    rownames(mat) <- NULL
    ht3 <- Heatmap(scale_rows(mat), show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression", cluster_columns=TRUE,show_column_names=FALSE,show_column_dend=FALSE,
             row_title="known",col = c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))) , 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE, #show_row_names=FALSE, 
            meth=anno_boxplot(mat, ylim=c(0,17), border=TRUE, outline=FALSE, pch=".", height= unit(1.5, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["known"], fontsize=6,gpar(fontsize=2)), width=unit(.2, "cm") )))

    #combine
    ht_list = ha %v% ht1 %v% ht2 %v% ht3
    #draw heatmap
    pdf(file.path(PostDE.dir, i, paste0("Complexheatmap_DEGup_stratClassCode_zScaleRaw.pdf")), height=12, width=8)
     draw(ht_list)
    dev.off()

    print(i)
}


col_fun = circlize::colorRamp2(quantile(scale_rows(assay(vst)[do.call("rbind",genes2plot$DACSB_vs_DMSO)$transcript_id ,]), c(seq(from=0, to=99, length=length(c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds")))))/100)),c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))))
for(i in names(genes2plot)){
    pheno <- colData(vst)[,]
    plot <- assay(vst)[genes2plot[[i]][[1]]$transcript_id,]
    ha <- HeatmapAnnotation(treatment = as.character(pheno$treatment),cellline= as.character(pheno$cellline), 
        col = list(treatment=treat_col[as.character(pheno$treatment)], cellline=cell_col[as.character(pheno$cellline)] ))

    #prepare heatmaps
    mat <- assay(vst)[genes2plot[[i]][["chimeric (novel)"]]$transcript_id,]
    rownames(mat) <- NULL
    ht1 <- Heatmap(scale_rows(mat), show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression", cluster_columns=TRUE,show_column_names=FALSE,show_column_dend=FALSE,
            row_title="chinmeric (novel)",col =col_fun, 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE,# show_row_names=FALSE, 
            meth=anno_boxplot(mat, ylim=c(0,17), border=TRUE, outline=FALSE, pch=".", height= unit(1.5, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["chimeric (novel)"],fontsize=6,annotation_nam=gpar(fontsize=6)), width=unit(.2, "cm") )))

    mat <- assay(vst)[genes2plot[[i]][["non-chimeric (novel)"]]$transcript_id,]
    rownames(mat) <- NULL
    ht2 <- Heatmap(scale_rows(mat), show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression",cluster_columns=TRUE,show_column_names=FALSE,show_column_dend=FALSE,
            row_title="non-chimeric (novel",col = col_fun, 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE,#show_row_names=FALSE, 
            meth=anno_boxplot(mat, ylim=c(0,17), border=TRUE, outline=FALSE, pch=".", height= unit(1.5, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["non-chimeric (novel)"], fontsize=6,annotation_nam=gpar(fontsize=6)), width=unit(.2, "cm") )))
    
     mat <- assay(vst)[genes2plot[[i]][["known"]]$transcript_id,]
    rownames(mat) <- NULL
    ht3 <- Heatmap(scale_rows(mat), show_row_dend=FALSE, name = "z Scaled\ntranscript\nexpression", cluster_columns=TRUE,show_column_names=FALSE,show_column_dend=FALSE,
             row_title="known",col = col_fun, 
            top_annotation=HeatmapAnnotation(show_annotation_name=FALSE, #show_row_names=FALSE, 
            meth=anno_boxplot(mat, ylim=c(0,17), border=TRUE, outline=FALSE, pch=".", height= unit(1.5, "cm"),gp = gpar(fill =treat_col[as.character(pheno$treatment)]))),
            left_annotation=rowAnnotation(cluster=anno_block( gpar(fill=class_col["known"], fontsize=6,gpar(fontsize=2)), width=unit(.2, "cm") )))

    #combine
    ht_list = ha %v% ht1 %v% ht2 %v% ht3
    #draw heatmap
    pdf(file.path(PostDE.dir, i, paste0("Complexheatmap_DEGup_stratClassCode_zScale.pdf")), height=12, width=8)
     draw(ht_list)
    dev.off()

    print(i)
}