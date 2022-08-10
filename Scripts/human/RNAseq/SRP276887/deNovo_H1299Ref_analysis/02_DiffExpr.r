#Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(knitr)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)
library(clusterProfiler)
library(readxl)
library(rtracklayer)

#function
plotMeta <- function(excel_path, n=10, output_path, color="#CD534CFF", height=3.5, width=3.5){
dataset <- read_excel(excel_path, sheet=2)
datasetsub <- dataset[grep("Summary", dataset$GroupID), ]
datasetsub <- head(datasetsub, n)
datasetsub$LogP<- abs(datasetsub$LogP)
datasetsub$Description<- rev(datasetsub$Description)
datasetsub$LogP<- rev(datasetsub$LogP)
datasetsub$Term<- rev(datasetsub$Term)
datasetsub$description <- as.factor(1:nrow(datasetsub)) 
 p <- ggplot(datasetsub, aes(x=description)) +
geom_bar(aes(x=description, y=LogP, fill = "-log10(p.value)"), stat = "identity", width=0.25)
 p<-p+ annotate("text",x=as.integer(datasetsub$description)+0.35, y=0.01, 
                label=datasetsub$Description, size=4, hjust = 0)
 p <- p + scale_y_continuous( expand = c(0, 0), name = "-log10(p-value)")
 p <- p + scale_fill_manual(values = c( color))
 p<-p+coord_flip() 
 p<-p+ ggpubr:::theme_pubr()
 p<-p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=15))
 p<-p+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +rremove(c("legend"))+rremove(c("ylab"))
 p$labels$fill <- ""
 pdf(output_path, height=height, width=width)
 print(p)
 dev.off()
}

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/220118_SRP276887_deNovo_quantification_H12Ref_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
dds <-readRDS(file.path(results.dir, "dds.rds"))
vst <- readRDS(file.path(results.dir, "vst.rds"))
dds <-readRDS(file.path(results.dir, "dds.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
#Take a look at design and annotation of samples
design(dds)
#Set specifications
alpha <- 0.05 #set FDR cutoff
lfc <- 1##set logfold2 cutoff

#Running the differential expression 
#set up all possible contrasts
contrasts <- as.data.frame(combn(as.character(unique(colData(dds)$treatment)), 2))
results <- list()
for (i in 1:length(contrasts)) {
  results[[i]]<- results(dds, contrast=c("treatment",as.character(contrasts[1,i]),  as.character(contrasts[2,i])), alpha = alpha, lfcThreshold = lfc) #extract results of wanted comparison
  print(paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i])))
  print(summary(results[[i]]))
}
nam<-vector()
for (i in 1:length(contrasts)) {
nam[i]<- paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i]))
}
names(results)<- nam
names(contrasts)<- nam


#annotate samples
#Make a dataframe out of it
DEG_results_list <- lapply(results, function(x){
  x<- as.data.frame(x)
  x$transcript_id <- rownames(x)
  x
})
names(DEG_results_list)

#add complete annotation
anno_transcripts <- anno[anno$type =="transcript",]
anno_genes_sub <- unique(as.data.frame(anno_transcripts)[, c("seqnames", "start", "end","strand","gene_id", "transcript_id", "dist_nearest_repeat", "nearest_repeat_repName", "nearest_repeat_repClass", "nearest_repeat_repFamily","nearest_LTR12repeat_repName","dist_nearest_LTR12repeat")])
DEG_results_list <- lapply(DEG_results_list, function(x){
  x <- dplyr::left_join(x, anno_genes_sub, by="transcript_id")
  rownames(x)<- x$transcript_id
  x
})

#add original reference annotation
anno_original_sub <- unique(as.data.frame(anno_original)[, c("transcript_id","transcript_type","gene_name" )])
DEG_results_list <- lapply(DEG_results_list, function(x){
  x <- dplyr::left_join(x, anno_original_sub, by="transcript_id")
  rownames(x)<- x$transcript_id
  x
})

#get average group levels
expre <- assay(vst)
pheno <- colData(vst)

DEG_results_list <- lapply(DEG_results_list, function(x){
  x$SETB1_KO <- NA
  x$SETB1_KO <- rowMeans(expre[, rownames(pheno[pheno$treatment=="SETB1_KO",])])
  x$Control <- NA
  x$Control <- rowMeans(expre[, rownames(pheno[pheno$treatment=="Control",])])
  x
})


#order by padjusted value
DEG_results_list <- lapply(DEG_results_list, function(x){
  x[which(is.na(x$padj)),]$padj <- 1
  x <- x[order(x$padj),]
  x
})
names(DEG_results_list) <- names(results)


#save lists
saveRDS(DEG_results_list, file.path(PostDE.dir, "DEG_results_group_list.rds"))
#DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#adjust daata for plotting
DEG_results_list_exp<- list()
for (i in names(DEG_results_list)){
    DEG_results_list_exp[[i]] <- as.data.frame(DEG_results_list[[i]])
    DEG_results_list_exp[[i]]$padj <- round(DEG_results_list_exp[[i]]$padj, 4)
    DEG_results_list_exp[[i]]$log2FoldChange <- round(DEG_results_list_exp[[i]]$log2FoldChange, 4)
    dir.create(file.path(PostDE.dir, i))
    write.table(DEG_results_list_exp[[i]], file=file.path(PostDE.dir,i, "DEG.tsv"),row.names =F, col.names = T , sep="\t", quote=FALSE)
    write.table(DEG_results_list_exp[[i]], file=file.path(PostDE.dir,i, "DEG.csv"),row.names =F, col.names = T , sep=", ", quote=FALSE)
}

#Plot results in Volcano and MA plot {.tabset}
DEG_results_list_plot <- DEG_results_list
DEG_results_list_plot <- lapply(DEG_results_list_plot, function(x){
  x <- as.data.frame(x)
  x$col <- "black"
  x$col <- ifelse(test = x$padj < alpha & x$log2FoldChange>lfc, yes = "red", 
                  no = ifelse(test = x$padj < alpha  & x$log2FoldChange<(-lfc), yes = "blue", "black"))
  temp <- c(head(x[x$col=="red",]$symbol, 5), head(x[x$col=="blue",]$symbol, 5))
  x$SYMBOL <- NA
  x[x$symbol %in% temp,]$SYMBOL <- x[x$symbol %in% temp,]$symbol
  x$pvalue <-(-(log10(x$pvalue)))
  x$padj <- (-(log10(x$padj)))
  x
}) 

lapply(DEG_results_list_plot, function(x){
    (table(x$col))
})
#Volcano Plot
library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- c("grey", col[3], col[1])

Volcanos<- list() 
for (i in names(DEG_results_list_plot)){
  Volcanos[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="log2FoldChange", y = "padj",
          xlab="Gene expression (log2 fold change)",#xscale="log10" ,
          ylab="- log10(P.adjusted)", #xlim=c(-0.5,0.5),
          color = "col" ,shape = 16, size = 1.5, palette=col, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance",# main=i# Customize reg. line
          #label="SYMBOL", repel = T, font.label = c( "black"), 
) + rremove("legend") + geom_vline(xintercept=c(-lfc,lfc), linetype = 2) +geom_hline(yintercept=-log10(alpha), linetype = 2)
  pdf(file.path(PostDE.dir,i,"volcano.pdf"))
  print(Volcanos[[i]])
  dev.off()
  pdf(file.path(PostDE.dir,i,"volcano2.pdf"), height=3.5, width=3.5)
  print(Volcanos[[i]])
  dev.off()
}

Volcanos<- list() 
for (i in names(DEG_results_list_plot)){
  Volcanos[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="log2FoldChange", y = "padj",
          xlab="Gene expression (log2 fold change)",#xscale="log10" ,
          ylab="- log10(P.adjusted)", #xlim=c(-0.5,0.5),
          color = "col" ,shape = 16, size = 1.5, palette=col, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance",# main=i# Customize reg. line
          label="SYMBOL", repel = T, font.label = c( "black"), 
) + rremove("legend") + geom_vline(xintercept=c(-lfc,lfc), linetype = 2) +geom_hline(yintercept=-log10(alpha), linetype = 2)
  pdf(file.path(PostDE.dir,i,"volcano3.pdf"), height=3.5, width=3.5)
  print(Volcanos[[i]])
  dev.off()
}

#MA plot
MAs <- list()
for (i in names(DEG_results_list_plot)){
  MAs[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="baseMean", y =  "log2FoldChange",xlab="Mean of normalized counts",xscale="log10" ,
          ylab="Gene expression (log2 fold change)",
          color = "col",palette=col,shape = 16, size = 1.5, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray") #legend.title="Significance",# Customize reg. line
          #label="SYMBOL", repel = TRUE#, legend = "right",
) +geom_hline(yintercept=c(-lfc,lfc), linetype = 2) + rremove("legend")
  pdf(file.path(PostDE.dir,i,"MA.pdf"))
  print(MAs[[i]])
  dev.off()
  pdf(file.path(PostDE.dir,i,"MA2.pdf"), height=3.5, width=3.5)
  print(MAs[[i]])
  dev.off()
}
MAs <- list()
for (i in names(DEG_results_list_plot)){
  MAs[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="baseMean", y =  "log2FoldChange",xlab="Mean of normalized counts",xscale="log10" ,
          ylab="Gene expression (log2 fold change)",
          color = "col",palette=col,shape = 16, size = 1.5, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), #legend.title="Significance",# Customize reg. line
          label="SYMBOL", repel = TRUE#, legend = "right",
) +geom_hline(yintercept=c(-lfc,lfc), linetype = 2) + rremove("legend")
  pdf(file.path(PostDE.dir,i,"MA3.pdf"), height=3.5, width=3.5)
  print(MAs[[i]])
  dev.off()
}

#heatmaps
#define padj cutoff for plotting
cutoff <- alpha
l2fc <- lfc
genes2plot <- lapply(DEG_results_list, function(x){
  x <- x[which(x$padj < cutoff & abs(x$log2FoldChange) >l2fc),]
  x <- rownames(x)
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]
lapply(genes2plot, function(x)length(x))

#heatmap
mat_genes<- assay(vst)
col <- c("#00AFBB", "#E7B800", "#FC4E07","gray")
names(col) <- levels(colData(dds)$treatment)
anno_colors <- list(treatment=col)
anno <- colData(vst)
annovst <- anno[,c("treatment" ), drop=FALSE]

heat<-list()
for (i in names(genes2plot)){
  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[[i]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
  annorow <- rownames(plot)
  heat[[i]]<- pheatmap(plot[, rownames(anno)], scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors
                       ) 
  pdf(file.path(PostDE.dir,i,"Heat_DEG_reds.pdf"))
  print(heat[[i]])
  dev.off()
    heat[[i]]<- pheatmap(plot[, rownames(anno)], scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors,
                       clustering_distance_rows="correlation") 
  pdf(file.path(PostDE.dir,i,"Heat_correlation_DEG_reds.pdf"))
  print(heat[[i]])
  dev.off()
}

#subset only groups of interest for heatmapa
heat<-list()
for (i in names(genes2plot)){
  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[[i]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
  annorow <- rownames(plot)

    comp_oi <- sapply(strsplit(i, "_vs_",fixed=TRUE), function(x)c(x[1], x[2]))
    heat[[i]]<- pheatmap(plot[,     rownames(anno[anno$treatment %in% comp_oi,])], 
    scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),show_rownames=F,annotation_colors=anno_colors,
                       clustering_distance_rows="correlation") 
  pdf(file.path(PostDE.dir,i,"Heat_correlation_DEG_subgroups.pdf"))
  print(heat[[i]])
  dev.off()
  print(i)
}

#common heatmapf
DEG <- unique(unlist(genes2plot))
length(DEG)
plot <- mat_genes[which(rownames(mat_genes) %in% DEG ),]
set.seed(123)
p <- pheatmap(plot, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
                       #cutree_rows =5, 
                       file=file.path(PostDE.dir, "Heat_allDEG_correlation.pdf")
                       ) 
p <- pheatmap(plot, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),
                       annotation_col=as.data.frame(annovst),
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
                       #cutree_rows =5, 
                       file=file.path(PostDE.dir, "Heat_allDEG_correlation_reds.pdf")
) 

# #metascape enrichment of the results of interest
# nrow(DEG_results_list$DACandSB939_vs_DMSO[which(DEG_results_list$DACandSB939_vs_DMSO$padj < alpha & DEG_results_list$DACandSB939_vs_DMSO$log2FoldChange>5),])
# paste0(DEG_results_list$DACandSB939_vs_DMSO[which(DEG_results_list$DACandSB939_vs_DMSO$padj < alpha & DEG_results_list$DACandSB939_vs_DMSO$log2FoldChange>5),]$gene_name, collapse=", ")


# #plot heatmap with clusters
# dir.create(file.path(PostDE.dir,"cluster_anal"))
# for(i in 4:8){
# #for(i in c(6)){
#     cut <- i
#     #cut tree
#     clust_genes <- sort(cutree(p$tree_row, k=cut))
#     clust_genes_df <- as.data.frame(clust_genes)-
#     colnames(clust_genes_df)<- "cluster"
#     clust_genes_df$cluster <- as.character(clust_genes_df$cluster)
#     #add color row annotation
#     color_cluster <- randomcoloR::distinctColorPalette(length(unique(clust_genes_df$cluster)))
#     names(color_cluster)<- unique(clust_genes_df$cluster)
#     anno_colors <- list(treatment=col, cluster=color_cluster)
#     #heatmap
#     p <- pheatmap(plot, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
#                         annotation_col=as.data.frame(annovst),
#                         annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
#                         cutree_rows =cut, annotation_row=clust_genes_df, 
#                         file=file.path(PostDE.dir,"cluster_anal", paste0("Heat_allDEG_correlation_clust", cut,".pdf")))
#     print(i)
# }

# # select 6 clusters
# cut <- 6
# #cut tree
# clust_genes <- sort(cutree(p$tree_row, k=cut))
# clust_genes_df <- as.data.frame(clust_genes)
# colnames(clust_genes_df)<- "cluster"
# clust_genes_df$cluster <- as.character(clust_genes_df$cluster)

# #write genes per cluster
# for (i in unique(clust_genes)){
# temp <- names(clust_genes[clust_genes==i])
# write.table(temp, file.path(PostDE.dir,"cluster_anal", paste0(cut,"cluster_clust",i, "_genes.txt")), quote=FALSE, row.names=FALSE)
# print(length(temp))
# }
# clust_genes<- as.data.frame(clust_genes)

# #plot results
# files <- list.files(file.path(PostDE.dir, "cluster_anal"), 
#     pattern="metascape_result.xlsx", full.names=TRUE, recursive=TRUE)
# paths <- gsub("metascape_result.xlsx", "", files)

# #plotting
# for(numb in c(3,5,10,15)){
# for(i in 1:length(files)){
#     plotMeta(files[i], 
#         file.path(paths[i], paste0("metascape_vis_cluster_",numb,"_1.pdf")), 
#         n=numb, width=5, height=5,color= "#86868699")
# }
# }

# #visualize average gene expression per cluster per group
# mat <- assay(vst)
# clust_expr_mean <- list()
# anno <- colData(vst)
# for(i in unique(clust_genes_df$cluster)){
#   genes_clust <- rownames(clust_genes_df[clust_genes_df$cluster==i,, drop=FALSE])
#   clust_expr<- mat[genes_clust,]
#   clust_expr_mean[[i]] <- colMeans(clust_expr)
#   clust_expr_mean[[i]]<- data.frame(AverageExpression=  clust_expr_mean[[i]], treatment=as.character(anno$treatment))
#   clust_expr_mean[[i]]$cluster <- paste0("cluster_", i)
# }
# clust_expr_mean_unlist <- do.call("rbind", clust_expr_mean)
# #order factors
# clust_expr_mean_unlist$cluster <- as.factor(clust_expr_mean_unlist$cluster)
# clust_expr_mean_unlist$cluster <- ordered(clust_expr_mean_unlist$cluster, paste0("cluster_", 1:6))

# pdf(file.path(PostDE.dir,"cluster_anal", "cluster_expression_boxplot.pdf"))
# ggboxplot(clust_expr_mean_unlist, x="treatment", y="AverageExpression", 
#   color="treatment", facet.by="cluster", palette=col) +
#   rotate_x_text(angle = 90) + stat_compare_means(label.y = .1, size=2.5)
# dev.off()
# #in 2 row visualization
# g <- ggboxplot(clust_expr_mean_unlist, x="treatment", y="AverageExpression", order=c("DMSO", "DAC", "SB939","DACandSB939"),
#   color="treatment", palette=col) +
#   rotate_x_text(angle = 90) + stat_compare_means(label.y = .1, size=2.5)
# g <- facet(g, facet.by="cluster", ncol=2)
# pdf(file.path(PostDE.dir,"cluster_anal", "cluster_expression_boxplot2.pdf"), height=14)
# g
# dev.off()
# g <- ggboxplot(clust_expr_mean_unlist, x="treatment", y="AverageExpression", order=c("DMSO", "DAC", "SB939","DACandSB939"),
#   color="treatment", palette=col) +
#   rotate_x_text(angle = 90)+ stat_compare_means(label.y = .1, size=2.5)
# g <- facet(g, facet.by="cluster", ncol=1)
# pdf(file.path(PostDE.dir,"cluster_anal", "cluster_expression_boxplot3_ordered.pdf"), height=20, width=7)
# g
# dev.off()

# #enrich for each cluster
# #bp
# GO<- list()
# for (i in unique(clust_genes$clust_genes)){
#     temp <- rownames(clust_genes[clust_genes==i,, drop=FALSE])
#     temp <- unique(sapply(strsplit(as.character(temp), ".", fixed=TRUE), "[",1))
#   GO[[i]] <- enrichGO(temp, OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
#                  pvalueCutoff = 1, qvalueCutoff = 1, ont = "BP")
#  print(i)
# }
# plot <- list()
# plot2 <- list()
# for (i in unique(clust_genes$clust_genes)){
#   plot[[i]] <- dotplot(GO[[i]],showCategory = 20)
#   pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_clust",i,"_GO_BP_dotplot.pdf")), width=12)
#   print(plot[[i]])
#   dev.off()
#   plot2[[i]] <- barplot(GO[[i]],showCategory = 10)
#   pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_clust",i,"_GO_BP_barplot.pdf")), height=4, width=10)
#   print(plot2[[i]])
#   dev.off()
# }
# #mf
# GO<- list()
# for (i in unique(clust_genes$clust_genes)){
#     temp <- rownames(clust_genes[clust_genes==i,, drop=FALSE])
#     temp <- unique(sapply(strsplit(as.character(temp), ".", fixed=TRUE), "[",1))
#   GO[[i]] <- enrichGO(temp, OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
#                  pvalueCutoff = 1, qvalueCutoff = 1, ont = "MF")
#  print(i)
# }
# plot <- list()
# plot2 <- list()
# for (i in unique(clust_genes$clust_genes)){
#   plot[[i]] <- dotplot(GO[[i]],showCategory = 20)
#   pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_clust",i,"_GO_MF_dotplot.pdf")), width=12)
#   print(plot[[i]])
#   dev.off()
#   plot2[[i]] <- barplot(GO[[i]],showCategory = 10)
#   pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_clust",i,"_GO_MF_barplot.pdf")), height=4, width=10)
#   print(plot2[[i]])
#   dev.off()
# }

# #kegg pathways
# kegg<- list()
# for (i in unique(clust_genes$clust_genes)){
#     temp <- rownames(clust_genes[clust_genes==i,, drop=FALSE])
#     temp <- unique(sapply(strsplit(as.character(temp), ".", fixed=TRUE), "[",1))
#     temp<- mapIds(org.Hs.eg.db, keys=temp, keytype ="ENSEMBL", column = "ENTREZID", multiVals = "first" )
#     kegg[[i]] <- enrichKEGG(temp, organism ="hsa",
#                     pvalueCutoff = 1, qvalueCutoff = 1)
#     print(i)
# }
# plot <- list()
# plot2 <- list()
# dir.create(file.path(PostDE.dir,"cluster_anal","enrich"))
# for (i in unique(clust_genes$clust_genes)){
#   plot[[i]] <- dotplot(kegg[[i]],showCategory = 20)
#   pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_clust",i,"_KEGG_dotplot.pdf")), width=12)
#   print(plot[[i]])
#   dev.off()
#   plot2[[i]] <- barplot(kegg[[i]],showCategory = 10)
#   pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_clust",i,"_KEGG_barplot.pdf")), height=4, width=10)
#   print(plot2[[i]])
#   dev.off()
# }




# #still to do

# #compare clusters
# clust_genes_df <- as.data.frame(clust_genes)
# colnames(clust_genes_df)<- "cluster"
# clust_genes_split <- split(clust_genes_df, clust_genes_df$cluster)
# clust_genes_split <- lapply(clust_genes_split, function(x)rownames(x))
# clust_genes_split<- lapply(clust_genes_split, function(x) mapIds(org.Hs.eg.db, keys=x, keytype ="SYMBOL", column = "ENTREZID", multiVals = "first" ))
# #Kegg
# formula_res <- compareCluster(clust_genes_split, fun="enrichKEGG", organism="mmu")#,pvalueCutoff = 1, qvalueCutoff = 1)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_KEGG_dotplot.pdf")), height=10, width=10)
# dotplot(formula_res) 
# dev.off()
# #GO MF
# formula_res <- compareCluster(clust_genes_split, fun="enrichGO", 
#     OrgDb = org.Hs.eg.db,keyType="ENTREZID", ont = "MF")
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_MF_dotplot.pdf")), height=10, width=11)
# dotplot(formula_res, includeAll=TRUE)
# dev.off()
# formula_res <- simplify(formula_res)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_MF_dotplot_simple.pdf")), height=10, width=11)
# dotplot(formula_res, includeAll=TRUE)
# dev.off()
# #GO BP
# formula_res <- compareCluster(clust_genes_split, fun="enrichGO", 
#     OrgDb = org.Hs.eg.db,keyType="ENTREZID", ont = "BP")
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_BP_dotplot.pdf")), height=10, width=11)
# dotplot(formula_res)
# dev.off()
# formula_res <- simplify(formula_res)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_BP_dotplot_simple.pdf")), height=10, width=11)
# dotplot(formula_res)
# dev.off()
# #GO CC
# formula_res <- compareCluster(clust_genes_split, fun="enrichGO", 
#     OrgDb = org.Hs.eg.db,keyType="ENTREZID", ont = "CC")
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_CC_dotplot.pdf")), height=10, width=12)
# dotplot(formula_res)
# dev.off()
# formula_res <- simplify(formula_res)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_CC_dotplot_simple.pdf")), height=10, width=12)
# dotplot(formula_res)
# dev.off()


# #same for stricter cutoffs
# #Kegg
# formula_res <- compareCluster(clust_genes_split, fun="enrichKEGG", organism="mmu", qvalueCutoff=0.01)#,pvalueCutoff = 1, qvalueCutoff = 1)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_KEGG_dotplot_strict.pdf")), height=10, width=10)
# dotplot(formula_res) 
# dev.off()
# #GO MF
# formula_res <- compareCluster(clust_genes_split, fun="enrichGO", 
#     OrgDb = org.Hs.eg.db,keyType="ENTREZID", ont = "MF", qvalueCutoff=0.01)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_MF_dotplot_strict.pdf")), height=10, width=11)
# dotplot(formula_res, includeAll=TRUE)
# dev.off()
# formula_res <- simplify(formula_res)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_MF_dotplot_simple_strict.pdf")), height=10, width=11)
# dotplot(formula_res, includeAll=TRUE)
# dev.off()
# #GO BP
# formula_res <- compareCluster(clust_genes_split, fun="enrichGO", 
#     OrgDb = org.Hs.eg.db,keyType="ENTREZID", ont = "BP", qvalueCutoff=0.01)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_BP_dotplot_strict.pdf")), height=10, width=11)
# dotplot(formula_res)
# dev.off()
# formula_res <- simplify(formula_res)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_BP_dotplot_simple_strict.pdf")), height=10, width=11)
# dotplot(formula_res)
# dev.off()
# #GO CC
# formula_res <- compareCluster(clust_genes_split, fun="enrichGO", 
#     OrgDb = org.Hs.eg.db,keyType="ENTREZID", ont = "CC", qvalueCutoff=0.01)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_CC_dotplot_strict.pdf")), height=10, width=12)
# dotplot(formula_res)
# dev.off()
# formula_res <- simplify(formula_res)
# pdf(file.path(PostDE.dir,"cluster_anal","enrich",paste0(cut,"cluster_Combined","_GO_CC_dotplot_simple_strict.pdf")), height=10, width=12)
# dotplot(formula_res)
# dev.off()

# #pies Distribution
# DEG_results_list_sub<- lapply(DEG_results_list, function(x){
#   x$direction <- NA
#   x$direction <- ifelse(x$log2FoldChange>0, "up", "down")
#   x
# })
# for(i in names(DEG_results_list_sub)){
#     temp <- DEG_results_list_sub[[i]][which(abs(DEG_results_list_sub[[i]]$log2FoldChange)>1 & DEG_results_list_sub[[i]]$padj <alpha),]
#     table <- as.data.frame(table(temp$direction))
#     labs <- paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq)

#     pdf(file.path(PostDE.dir,i,paste0("Direction_DEGs_sign_fc.pdf")), height=3.5, width=3.5)
#     print(ggpie(table, "Freq", fill="Var1", palette=col[c(2,3)], label=labs, lab.pos = "in", main="DEG Analysis\nadj. p-value < alpha &\nabs. log2 fold-change > 1",  lab.font = c(5, "bold", "white")) + rremove("legend"))
#     dev.off()

#     temp <- DEG_results_list_sub[[i]][which(DEG_results_list_sub[[i]]$padj <alpha),]
#     table <- as.data.frame(table(temp$direction))
#     labs <- paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq)

#     pdf(file.path(PostDE.dir,i,paste0("Direction_DEGs_sign.pdf")), height=3.5, width=3.5)
#     print(ggpie(table, "Freq", fill="Var1", palette=col[c(2,3)], label=labs, lab.pos = "in", main="DEG Analysis\nadj. p-value < alpha",  lab.font = c(5, "bold", "white")) + rremove("legend"))
#     dev.off()
# }




