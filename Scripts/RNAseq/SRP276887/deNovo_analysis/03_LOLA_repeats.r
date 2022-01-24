#Perform differential epression analysis of different contrasts
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
#function 
label_func <- function(x){
    breaks <- x
    #breaks[breaks>=500] <-  ">=500"
    #breaks <- droplevels(breaks)
    breaks

}
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare_SRP276887_mergedTranscripts.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
#Take a look at design and annotation of samples
design(dds)
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 1##set logfold2 cutoff

#load database
regionDB_repFamily_hg19 <- loadRegionDB(file.path(datasets.dir,"repFamily_hg19"))
regionDB_repClass_hg19 <- loadRegionDB(file.path(datasets.dir,"repClass_hg19"))
regionDB_repName_LTR_hg19 <- loadRegionDB(file.path(datasets.dir,"repName_LTR_hg19"))


#run enrichment analysis
anno_transcript <- anno[anno$type == "transcript",]
userUnisverse <- resize(anno_transcript,1)
results_repFamily_hg19 <- list()
results_repClass_hg19 <- list()
results_repName_LTR_hg19 <-list()
for(i in names(DEG_results_list)){
    #get start position(strand specific)
    UserSet <- list(Upregulated=resize(anno_transcript[anno_transcript$transcript_id %in%rownames(DEG_results_list[[i]][which(DEG_results_list[[i]]$padj < alpha & DEG_results_list[[i]]$log2FoldChange>lfc),]),],1),
        Downregulataed=resize(anno_transcript[anno_transcript$transcript_id %in%rownames(DEG_results_list[[i]][which(DEG_results_list[[i]]$padj < alpha & DEG_results_list[[i]]$log2FoldChange<(-lfc)),]),],1))
    results_repFamily_hg19[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repFamily_hg19, cores=4)
    results_repClass_hg19[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repClass_hg19, cores=4)
    results_repName_LTR_hg19[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repName_LTR_hg19, cores=4)
    print(i)
}
dir.create(file.path(PostDE.dir, "LOLA"))
saveRDS(results_repFamily_hg19, file.path(PostDE.dir, "LOLA", "results_repFamily_hg19.rds"))
saveRDS(results_repClass_hg19, file.path(PostDE.dir, "LOLA", "results_repClass_hg19.rds"))
saveRDS(results_repName_LTR_hg19, file.path(PostDE.dir, "LOLA", "results_repName_LTR_hg19.rds"))
#results_repFamily_hg19 <- readRDS(file.path(PostDE.dir, "LOLA", "results_repFamily_hg19.rds"))
#results_repClass_hg19 <- readRDS(file.path(PostDE.dir, "LOLA", "results_repClass_hg19.rds"))
#results_repName_LTR_hg19 <- readRDS(file.path(PostDE.dir, "LOLA", "results_repName_LTR_hg19.rds"))


#Plot genomic regions in bubble plot
#subset result of interest
for(i in names(DEG_results_list)){
dir.create(file.path(PostDE.dir, i,"LOLA"))

results <- results_repFamily_hg19[[i]]
#plotting
idx <- head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i,"LOLA", "results_repFamily_LTR_hg19_top20.pdf"), height=10)
print(g)
dev.off()
#all
#plotting
idx <- head(unique(results$filename),Inf)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA", "results_repFamily_LTR_hg19_topInf.pdf"), height=20)
print(g)
dev.off()

#plot sign ones
idx <- unique(results[results$qValue < 0.05, ]$filename)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA", "results_repFamily_LTR_hg19_Sign.pdf"), height=20)
print(g)
dev.off()


#Plot genomic regions in bubble plot
#subset result of interest
results <- results_repClass_hg19[[i]]
#plotting
idx <- head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA", "results_repClass_LTR_hg19_top20.pdf"), height=10)
print(g)
dev.off()
#all
#plotting
idx <- head(unique(results$filename),Inf)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA", "results_repClass_LTR_hg19_topInf.pdf"), height=20)
print(g)
dev.off()

#plot sign ones
idx <- unique(results[results$qValue < 0.05, ]$filename)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA", "results_repClass_LTR_hg19_Sign.pdf"), height=20)
print(g)
dev.off()

#Plot genomic regions in bubble plot
#subset result of interest
results <- results_repName_LTR_hg19[[i]]
#plotting
idx <- head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA", "results_repName_LTR_hg19_top20.pdf"), height=10)
print(g)
dev.off()
idx <- head(unique(results$filename),5)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA","results_repName_LTR_hg19_top5.pdf"), height=10)
print(g)
dev.off()
#all
#plotting
idx <- head(unique(results$filename),Inf)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA","results_repName_LTR_hg19_topInf.pdf"), height=20)
print(g)
dev.off()

#plot sign ones
idx <- unique(results[results$qValue < 0.05, ]$filename)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="blue", high="red", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path(PostDE.dir, i, "LOLA","results_repName_LTR_hg19_Sign.pdf"), height=20)
print(g)
dev.off()

print(i)
}



