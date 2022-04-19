#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)

#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Extracting transformed values
vst <- readRDS(file.path(results.dir,"vst.rds"))
anno <- colData(vst)

#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c( "DMSO", "SB939","DAC", "DACSB")
comp_col <- col[c(2,8,6)]
names(comp_col)<-c( "DAC_vs_DMSO", "SB939_vs_DMSO",  "DACandSB939_vs_DMSO")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")
ltr_col <- brewer.pal(5, "Accent")
names(ltr_col)<- c("LTR12","LTR12C","LTR12D","LTR12F", "no LTR12")
cell_col  <- brewer.pal(length(unique(colData(vst)$cellline)), "Paired")
names(cell_col)<- as.character(unique(colData(vst)$cellline))

#run pca
set.seed(43)
topVar<- head(order(rowVars(assay(vst)), decreasing=TRUE),  5000)
matvst <- assay(vst)[topVar,]

ir.pca <- prcomp(t(as.matrix(matvst)),
                 center = T,
                 scale. = F) 
    
print(summary(ir.pca))
pc <- ir.pca$x
pc<- as.data.frame(cbind(pc ,anno))
pc$treatment <- as.character(pc$treatment)

dir.create(file.path(PreDE.dir, "PCA"))
pdf(file.path(PreDE.dir, "PCA","PCA12_5000mvgenes.pdf"), height = 5, width = 6)
ggscatter(pc, x="PC1", y="PC2",
            color = "treatment", 
            shape= "cellline",
            ellipse = T,, mean.point = FALSE,palette=treat_col,
            star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), "% variance")), 
            ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance")))+
            theme(legend.position="right")
dev.off()
pdf(file.path(PreDE.dir, "PCA","PCA13_5000mvgenes.pdf"), height = 5, width = 6)
ggscatter(pc, x="PC1", y="PC3",
            color = "treatment", 
            shape= "cellline",
            ellipse = F, mean.point = FALSE,palette=treat_col,
            star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), "% variance")), 
            ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance")))+
            theme(legend.position="right")
dev.off()
pdf(file.path(PreDE.dir, "PCA","PCA13_5000mvgenes_withEllipse.pdf"), height = 5, width = 6)
ggscatter(pc, x="PC1", y="PC3",
            color = "treatment", 
            shape= "cellline",
            ellipse = T,, mean.point = FALSE,palette=treat_col,
            star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), "% variance")), 
            ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance")))+
            theme(legend.position="right")
dev.off()