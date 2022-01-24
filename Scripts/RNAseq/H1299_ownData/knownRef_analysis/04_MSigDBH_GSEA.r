#perform gsea with msigdb hallmarks

#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(knitr)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)
library(clusterProfiler)
library(readxl)
library(rtracklayer)
library(LOLA)
library(ChIPpeakAnno)
library(rtracklayer)
library(ggpubr)
library(LOLA)
library(GenomicFeatures)
library(clusterProfiler)

#function 
plotGSEA <- function(input, n, path, height=7, width=7){
    SainiGSE.df  <- as.data.frame(input)
    SainiGSE.df <- SainiGSE.df [order(abs(SainiGSE.df$NES), decreasing = TRUE),]
    SainiGSEsub <- SainiGSE.df[, c("NES", "p.adjust","ID", "pvalue")]
    SainiGSEsub<- head(SainiGSEsub, n)
    SainiGSEsub$pvalue <- abs(log10(SainiGSEsub$p.adjust ))
    SainiGSEsub$factor <- as.factor(c(1:nrow(SainiGSEsub)))
    SainiGSEsub$factor<- rev(SainiGSEsub$factor)
    p<-ggplot(SainiGSEsub, aes(x=factor)) +
    geom_bar(aes(x=as.integer(as.factor(SainiGSEsub$factor))-0.15, y=NES, fill = "Normalized Enrichment Score"), stat = "identity", width=0.25)+ 
    geom_bar(aes(x=as.integer(as.factor(SainiGSEsub$factor))+0.15, y=pvalue, fill = "p.value"), stat = "identity", width=0.25)
    SainiGSEsub$ID <- gsub("."," ",SainiGSEsub$ID, fixed=TRUE)
    SainiGSEsub$ID <- gsub("  "," ",SainiGSEsub$ID, fixed=TRUE)
    SainiGSEsub$ID <- gsub("_"," ",SainiGSEsub$ID, fixed=TRUE)
    SainiGSEsub$ID <- gsub("HALLMARK","",SainiGSEsub$ID, fixed=TRUE)


    p<-p+ annotate("text",x=as.integer(as.factor(SainiGSEsub$factor))+0.5, y=0.01, 
                label=SainiGSEsub$ID, size=5, hjust = 0)
    p <- p + scale_y_continuous(sec.axis = sec_axis(~.*2, name = "-log10(adj. p-value)"), expand = c(0, 0))
    p <- p + scale_fill_manual(values = c("#CD534CFF", "grey"))
    p<-p+coord_flip() 
    p<-p+ ggpubr:::theme_pubr()
    p<-p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=15))
    p<-p+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
    p <- p + labs(y = "Normalized enrichment score",
                x = NULL,
                                colour = "Parameter")
    p$labels$fill <- ""
    pdf(path, height=height, width=width)
    print(p)
    dev.off()
}
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
ext.dir <- "/omics/groups/OE0219/internal/tinat/external_data/"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
anno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")

#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

#Import hallmark Gene Sets
H <- read.gmt(file.path(ext.dir, "misgdb", "h.all.v7.5.symbols.gmt"))

#generate list of degs
genelists <- lapply(DEG_results_list, function(x){
    y <- x[, "log2FoldChange"]
    names(y)<- x$gene_name
    y <- sort(y, decreasing=TRUE)
    y
})

#run gsea
gsea_H<-list()
for(i in names(genelists)){
    gsea_H[[i]] <- GSEA(genelists[[i]], TERM2GENE=H, pvalueCutoff=1)
    dir.create(file.path(PostDE.dir, i, "MSigDB"))
    plotGSEA(input=gsea_H[[i]], n=10, path=file.path(PostDE.dir, i, "MSigDB","MSigDB_Hallmarks_GSEA_results_barplot.pdf"), height=7, width=7)

}


#Write GSE barcode plots
for ( i in c("DACandSB939_vs_DMSO", "DAC_vs_DMSO","SB939_vs_DMSO")){
    temp <- gsea_H[[i]]
    for(j in 1:length(temp$Description)){
            pdf(file.path(PostDE.dir, i,  "MSigDB",paste0("MSigDB_Hallmarks_GSEA_",temp$Description[[j]], "_barcode.pdf")),width = 7,height=5)
                print(gseaplot(temp, geneSetID =temp$ID[j] ,color.line ="#CD534CFF",
                    title=paste0("NES=",round(temp$NES[j],3),"; adj. p-value=",round(temp$p.adjust[j],3))))
            dev.off()
    }
    print(i)
}



#same for CTA
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

#Import hallmark Gene Sets
CTA <- read.gmt(file.path(ext.dir, "misgdb", "CTA_geneset.gmt"))

#generate list of degs
genelists <- lapply(DEG_results_list, function(x){
    y <- x[, "log2FoldChange"]
    names(y)<- x$gene_name
    y <- sort(y, decreasing=TRUE)
    y
})

#run gsea
gsea_H<-list()
for(i in names(genelists)){
    gsea_H[[i]] <- GSEA(genelists[[i]], TERM2GENE=CTA, pvalueCutoff=1)
    #plotGSEA(input=gsea_H[[i]], n=10, path=file.path(PostDE.dir, i, "MSigDB","MSigDB_Hallmarks_GSEA_results_barplot.pdf"), height=7, width=7)

}


#Write GSE barcode plots
for ( i in c("DACandSB939_vs_DMSO", "DAC_vs_DMSO","SB939_vs_DMSO")){
    temp <- gsea_H[[i]]
    for(j in 1:length(temp$Description)){
    #for(j in c("HALLMARK_INTERFERON_ALPHA_RESPONSE")){
            pdf(file.path(PostDE.dir, i,  "MSigDB",paste0("MSigDB_Hallmarks_GSEA_",temp$Description[[j]], "_barcode.pdf")),width = 7,height=5)
                print(gseaplot(temp, geneSetID =temp$ID[j] ,color.line ="#CD534CFF",
                    title=paste0(temp$ID[j], "\nNES=",round(temp$NES[j],3),"; adj. p-value=",round(temp$p.adjust[j],3))))
            dev.off()

    }
    print(i)
}

#boxplto of testis gene expression
anno_sub <- anno_original[anno_original$gene_name %in% CTA[,2]]
mat <- assay(vst)
mat_sub <- mat[rownames(mat) %in% anno_sub$gene_id,]
mat_sub_treatment <- data.frame(DMSO=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="DMSO",])]), 
    DAC=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="DAC",])]),
    SB939=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="SB939",])]),
    DACandSB939=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="DACandSB939",])]))

mat_sub_melt <- reshape2::melt(mat_sub_treatment)
col <- c("#00AFBB", "#E7B800", "#FC4E07","gray")
names(col) <- levels(colData(dds)$treatment)
pdf(file.path(PostDE.dir,paste0("CTA_expression_levels.pdf")),width = 5,height=5)
ggboxplot(mat_sub_melt, x="variable", y="value", ylab="Average CTA expression", fill="variable", palette=col)+
    rremove("xlab")+rremove("legend")
dev.off()

#same for aim
AIM <- as.data.frame(data.table::fread(file.path(ext.dir, "Li_etal", "ubiquitous_AIMS.txt")))
CTA <-  as.data.frame(data.table::fread(file.path(ext.dir,  "Cancer_Testis_Antigens.txt"), header=FALSE))
AIM$ont <- "AZA IMmune gene set"
CTA$ont <- "Cancer testis antigens"
GS <- data.frame(ont=c(AIM$ont, CTA$ont), gene=c(AIM$ADM, CTA$V1))


#run gsea
gsea_H<-list()
for(i in names(genelists)){
    gsea_H[[i]] <- GSEA(genelists[[i]], TERM2GENE=GS, pvalueCutoff=1)
    plotGSEA(input=gsea_H[[i]], n=10, path=file.path(PostDE.dir, i,"AIM_CTA_GSEA_results_barplot.pdf"), height=7, width=7)

}

#Write GSE barcode plots
for ( i in c("DACandSB939_vs_DMSO", "DAC_vs_DMSO","SB939_vs_DMSO")){
    temp <- gsea_H[[i]]
    for(j in 1:length(temp$Description)){
            pdf(file.path(PostDE.dir, i,  "MSigDB",paste0("MSigDB_Hallmarks_GSEA_",temp$Description[[j]], "_barcode.pdf")),width = 7,height=5)
                print(gseaplot(temp, geneSetID =temp$ID[j] ,color.line ="#CD534CFF",
                    title=paste0(temp$ID[j], "\nNES=",round(temp$NES[j],3),"; adj. p-value=",round(temp$p.adjust[j],3))))
            dev.off()

    }
    print(i)
}

#boxplto of testis gene expression
anno_sub <- anno_original[anno_original$gene_name %in% AIM[,2]]
mat <- assay(vst)
mat_sub <- mat[rownames(mat) %in% anno_sub$gene_id,]
mat_sub_treatment <- data.frame(DMSO=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="DMSO",])]), 
    DAC=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="DAC",])]),
    SB939=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="SB939",])]),
    DACandSB939=rowMeans(mat_sub[,rownames(pheno[pheno$treatment=="DACandSB939",])]))

mat_sub_melt <- reshape2::melt(mat_sub_treatment)
col <- c("#00AFBB", "#E7B800", "#FC4E07","gray")
names(col) <- levels(colData(dds)$treatment)
pdf(file.path(PostDE.dir,paste0("AIM_expression_levels.pdf")),width = 5,height=5)
ggboxplot(mat_sub_melt, x="variable", y="value", ylab="Average AIM expression", fill="variable", palette=col)+
    rremove("xlab")+rremove("legend")
dev.off()