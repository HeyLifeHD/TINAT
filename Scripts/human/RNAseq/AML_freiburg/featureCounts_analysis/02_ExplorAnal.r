#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)
#function
plotPCA.alt <- function (object, intgroup = "condition", ntop = Inf, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + 
    geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220622_AMLPatient_FeatureCountRepeat_analysis_revised/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds_repFamily_group.rds"))

#Extracting transformed values
vst <- vst(dds)
anno <- colData(vst)
saveRDS(vst,file = file.path(results.dir,"vst.rds"))
#vst <- readRDS(file.path(results.dir,"vst.rds"))


#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2)]
#names(treat_col)<-c( "DMSO", "SB939","DAC", "DACSB")
names(treat_col)<-c( "pre_treatment", "post_treatment")
#comp_col <- col[c(2,8,6)]
#names(comp_col)<-c( "DAC_vs_DMSO", "SB939_vs_DMSO",  "DACandSB939_vs_DMSO")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
ltr_col <- brewer.pal(5, "Accent")
names(ltr_col)<- c("LTR12","LTR12C","LTR12D","LTR12F", "no LTR12")
#cell_col  <- brewer.pal(length(unique(colData(vst)$Patient_ID)), "Paired")
#names(cell_col)<- as.character(unique(colData(vst)$Patient_ID))
cell_col  <- c(brewer.pal(10, "Paired"), brewer.pal(4, "Dark2"))
names(cell_col)<- as.character(unique(colData(vst)$Patient_ID))
#anno
anno_colors <- list(group=treat_col, Patient_ID= cell_col)
anno <- colData(vst)
annovst <- anno[,c("Patient_ID", "group" ), drop=FALSE]

#plot pca
pcaDatavst<- plotPCA(vst, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarvst <- round(100* attr(pcaDatavst, "percentVar"))
names(col) <- levels(colData(dds)$group)
pdf(file.path(PreDE.dir, "PCA12.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst, x="PC1", y="PC2",
            size=5,
            color = "group", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst[2], "% variance")) )
dev.off()


#pca of 5000 mv genes
pcaDatavst<- plotPCA(vst, intgroup =  colnames(anno), returnData =TRUE, ntop=5000)
percentvarvst <- round(100* attr(pcaDatavst, "percentVar"))
pdf(file.path(PreDE.dir, "PCA12_5000mvGenes.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst, x="PC1", y="PC2",
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst[2], "% variance")) )
dev.off()
pcaDatavst<- plotPCA(vst, intgroup =  colnames(anno), returnData =TRUE, ntop=1000)
percentvarvst <- round(100* attr(pcaDatavst, "percentVar"))
pdf(file.path(PreDE.dir, "PCA12_1000mvGenes.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst, x="PC1", y="PC2",
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst[2], "% variance")) )
dev.off()

#function to plot PC2 and 3
pcaDatavst.alt<- plotPCA.alt(vst, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
percentvarvst.alt <- round(100* attr(pcaDatavst.alt, "percentVar"))
pdf(file.path(PreDE.dir, "PCA23.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst.alt, x="PC2", y="PC3",            
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
          ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarvst.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarvst.alt[2], "% variance")) )
dev.off()

#pca 2 and 3 of 5000 mv genes
pcaDatavst.alt<- plotPCA.alt(vst, intgroup = colnames(anno), returnData =TRUE, ntop=5000)
percentvarvst.alt <- round(100* attr(pcaDatavst.alt, "percentVar"))
pdf(file.path(PreDE.dir, "PCA23_5000mvGenes.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst.alt, x="PC2", y="PC3",            
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
          ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarvst.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarvst.alt[2], "% variance")) )
dev.off()

#Sample Clustering correlation
mat <- assay(vst)
dissimilarity <- 1 - cor(mat, use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
col1 <- hrld$labels
names(col1)<- as.character(anno[col1,]$group)
col1 <- col1[order.dendrogram(dend)]
col1 <- treat_col[names(col1)]
labels(dend)<- as.character(anno[labels(dend),]$Patient_ID)
dend <- dend %>% 
set("branches_lwd", 2) %>%
set("labels_colors",col1) %>% 
set("labels_cex", .6 )%>%
set("leaves_pch", 19)%>% 
set("leaves_cex", 1.5)%>% 
set("leaves_col", col1)
pdf(file.path(PreDE.dir, "Clustering_correlation.pdf"), height = 5, width = 5)
dend %>% plot
dev.off()

#Gene clustering: Heatmaps
#500 most variable repeats
topVarrepeatsvst<- head(order(rowVars(assay(vst)), decreasing=TRUE),  1000)
matvst <- assay(vst)[topVarrepeatsvst,]
pdf(file.path(PreDE.dir,"Heatmap1000vst_Scale_repeats.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors)
dev.off()

topVarrepeatsvst<- head(order(rowVars(assay(vst)), decreasing=TRUE),  100)
matvst <- assay(vst)[topVarrepeatsvst,]
pdf(file.path(PreDE.dir,"Heatmap100vst_Scale_repeats.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

#only ltr
vst_LTR <-  assay(vst)[grep("LTR12",rownames( assay(vst))),]

pdf(file.path(PreDE.dir,"Heatmap_noScale_LTR12repeats.pdf"),height= 7)
pheatmap(vst_LTR, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(vst_LTR),
annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors)
dev.off()


vst_LTR <-  assay(vst)[rownames( assay(vst)) %in% c("LTR12C", "LTR12D"),]

pdf(file.path(PreDE.dir,"Heatmap_noScale_LTR12CDrepeats.pdf"),height= 7)
pheatmap(vst_LTR, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(vst_LTR),
annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors)
dev.off()







#redo with removed batch effect
modmat <- model.matrix(~ Patient_ID + group ,as.data.frame(colData(vst))[, c("Patient_ID", "group")])

BR<-removeBatchEffect(assay(vst), batch=as.vector(colData(vst)$Patient_ID), design=modmat)
vst_b <- SummarizedExperiment(assays =BR , colData=colData(vst))
vst_b <- DESeqTransform(vst_b)

saveRDS(vst_b,file =file.path(results.dir, "vst_b.rds"))
vst_b <- readRDS(file.path(results.dir, "vst_b.rds"))
 

#plot pca
pcaDatavst_b<- plotPCA(vst_b, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarvst_b <- round(100* attr(pcaDatavst_b, "percentVar"))
col <- c("#00AFBB", "#E7B800", "#FC4E07", "gray")
names(col) <- levels(colData(dds)$group)
pdf(file.path(PreDE.dir, "PCA12_batchRemoved.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst_b, x="PC1", y="PC2",
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst_b[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst_b[2], "% variance")) )
dev.off()

#pca of 5000 mv genes
pcaDatavst_b<- plotPCA(vst_b, intgroup =  colnames(anno), returnData =TRUE, ntop=5000)
percentvarvst_b <- round(100* attr(pcaDatavst_b, "percentVar"))
pdf(file.path(PreDE.dir, "PCA12_5000mvGenes_batchRemoved.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst_b, x="PC1", y="PC2",
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst_b[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst_b[2], "% variance")) )
dev.off()
pcaDatavst_b<- plotPCA(vst_b, intgroup =  colnames(anno), returnData =TRUE, ntop=1000)
percentvarvst_b <- round(100* attr(pcaDatavst_b, "percentVar"))
pdf(file.path(PreDE.dir, "PCA12_1000mvGenes_batchRemoved.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst_b, x="PC1", y="PC2",
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarvst_b[1], "% variance")), ylab=(paste0("PC2: ", percentvarvst_b[2], "% variance")) )
dev.off()

#function to plot PC2 and 3
pcaDatavst_b.alt<- plotPCA.alt(vst_b, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
percentvarvst_b.alt <- round(100* attr(pcaDatavst_b.alt, "percentVar"))
pdf(file.path(PreDE.dir, "PCA23_batchRemoved.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst_b.alt, x="PC2", y="PC3",            
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
          ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarvst_b.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarvst_b.alt[2], "% variance")) )
dev.off()

#pca 2 and 3 of 5000 mv genes
pcaDatavst_b.alt<- plotPCA.alt(vst_b, intgroup = colnames(anno), returnData =TRUE, ntop=5000)
percentvarvst_b.alt <- round(100* attr(pcaDatavst_b.alt, "percentVar"))
pdf(file.path(PreDE.dir, "PCA23_5000mvGenes_batchRemoved.pdf"), height = 5, width = 6)
ggscatter(pcaDatavst_b.alt, x="PC2", y="PC3",            
            size=5,
            color = "group.1", 
            shape= "Patient_ID",
            #label= "replicate",
            repel=TRUE,
            legend = "right",palette=treat_col,
          ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarvst_b.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarvst_b.alt[2], "% variance")) )
dev.off()

#Sample Clustering correlation
mat <- assay(vst_b)
dissimilarity <- 1 - cor(mat, use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
col1 <- hrld$labels
names(col1)<- as.character(anno[col1,]$group)
col1 <- col1[order.dendrogram(dend)]
col1 <- treat_col[names(col1)]
labels(dend)<- as.character(anno[labels(dend),]$Patient_ID)
dend <- dend %>% 
set("branches_lwd", 2) %>%
set("labels_colors",col1) %>% 
set("labels_cex", .6 )%>%
set("leaves_pch", 19)%>% 
set("leaves_cex", 1.5)%>% 
set("leaves_col", col1)
pdf(file.path(PreDE.dir, "Clustering_correlation_batchRemoved.pdf"), height = 5, width = 5)
dend %>% plot
dev.off()

#Sample Clustering correlation
topVar<- head(order(rowVars(assay(vst_b)), decreasing=TRUE),  1000)
matvst_b <- assay(vst_b)[topVar,]
dissimilarity <- 1 - cor(matvst_b, use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)
hrld<- hclust(distance)
dend<- hrld%>% as.dendrogram 
col1 <- hrld$labels
names(col1)<- as.character(anno[col1,]$group)
col1 <- col1[order.dendrogram(dend)]
col1 <- treat_col[names(col1)]
labels(dend)<- as.character(anno[labels(dend),]$Patient_ID)
dend <- dend %>% 
set("branches_lwd", 2) %>%
set("labels_colors",col1) %>% 
set("labels_cex", .6 )%>%
set("leaves_pch", 19)%>% 
set("leaves_cex", 1.5)%>% 
set("leaves_col", col1)
pdf(file.path(PreDE.dir, "Clustering_correlation_1000mvTrans_batchRemoved.pdf"), height = 5, width = 5)
dend %>% plot
dev.off()

#Gene clustering: Heatmaps
#500 most variable repeats
topVarrepeatsvst<- head(order(rowVars(assay(vst_b)), decreasing=TRUE),  1000)
matvst <- assay(vst)[topVarrepeatsvst,]
pdf(file.path(PreDE.dir,"Heatmap1000vst_Scale_repeats_batchRemoved.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors)
dev.off()

topVarrepeatsvst<- head(order(rowVars(assay(vst_b)), decreasing=TRUE),  100)
matvst <- assay(vst)[topVarrepeatsvst,]
pdf(file.path(PreDE.dir,"Heatmap100vst_Scale_repeats_batchRemoved.pdf"),height= 7)
pheatmap(matvst, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(matvst),
                       annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors, fontsize_row=5)
dev.off()

#only ltr
vst_LTR <-  assay(vst)[grep("LTR12",rownames( assay(vst_b))),]

pdf(file.path(PreDE.dir,"Heatmap_noScale_LTR12repeats_batchRemoved.pdf"),height= 7)
pheatmap(vst_LTR, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(vst_LTR),
annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors)
dev.off()


vst_LTR <-  assay(vst)[rownames( assay(vst_b)) %in% c("LTR12C", "LTR12D"),]

pdf(file.path(PreDE.dir,"Heatmap_noScale_LTR12CDrepeats_batchRemoved.pdf"),height= 7)
pheatmap(vst_LTR, scale="row", show_colnames=F,color= c(hcl.colors(20,"Blues"),rev(hcl.colors(20,"Reds"))),labels_row=rownames(vst_LTR),
annotation_col=as.data.frame(annovst), show_rownames=T,annotation_colors=anno_colors)
dev.off()
