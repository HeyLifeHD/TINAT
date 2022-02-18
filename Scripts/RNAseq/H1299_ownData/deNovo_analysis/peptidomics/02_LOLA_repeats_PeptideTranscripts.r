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
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
anno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
#Take a look at design and annotation of samples
design(dds)
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

#load database
regionDB_repFamily_hg19 <- loadRegionDB(file.path(datasets.dir,"repFamily_hg19"))
regionDB_repClass_hg19 <- loadRegionDB(file.path(datasets.dir,"repClass_hg19"))
regionDB_repName_LTR_hg19 <- loadRegionDB(file.path(datasets.dir,"repName_LTR_hg19"))

#load peptide data
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))

#subset peptide list based on our ORFs == Universe
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]

#get transcripts that give rise to uniquely presented peptides
peptides_new_ORF_noDMSO <- peptides_new_ORF[which(peptides_new_ORF$DMSO ==0),]
peptides_new_ORF_noDMSO_split <- split(peptides_new_ORF_noDMSO, peptides_new_ORF_noDMSO$DAC_SB)
lapply(peptides_new_ORF_noDMSO_split, length)
names(peptides_new_ORF_noDMSO_split)<- c("1/3 DAC + SB939 unique", "2/3 DAC + SB939 unique", "3/3 DAC + SB939 unique")

#only get transcript names 
combined_analysis <- list(combined_analysis=lapply(peptides_new_ORF_noDMSO_split, function(x){
    x <- unique(x$transcript_id)
    x <- resize(anno_transcript[anno_transcript$transcript_id %in%x,],1)
    x
    }))


#run enrichment analysis
anno_transcript <- anno[anno$type == "transcript",]
userUnisverse <- resize(anno_transcript[anno_transcript$transcript_id %in% peptides_new_ORF$transcript_id],1)
results_repFamily_hg19 <- list()agne
results_repClass_hg19 <- list()
results_repName_LTR_hg19 <-list()
for(i in names(combined_analysis)){
    #get start position(strand specific)
    UserSet <- combined_analysis[[i]]
    results_repFamily_hg19[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repFamily_hg19, cores=4)
    results_repClass_hg19[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repClass_hg19, cores=4)
    results_repName_LTR_hg19[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repName_LTR_hg19, cores=4)
    print(i)
}
dir.create(file.path(output.dir, "LOLA"))
saveRDS(results_repFamily_hg19, file.path(output.dir, "LOLA", "results_repFamily_hg19.rds"))
saveRDS(results_repClass_hg19, file.path(output.dir, "LOLA", "results_repClass_hg19.rds"))
saveRDS(results_repName_LTR_hg19, file.path(output.dir, "LOLA", "results_repName_LTR_hg19.rds"))
#results_repFamily_hg19 <- readRDS(file.path(output.dir, "LOLA", "results_repFamily_hg19.rds"))
#results_repClass_hg19 <- readRDS(file.path(output.dir, "LOLA", "results_repClass_hg19.rds"))
#results_repName_LTR_hg19 <- readRDS(file.path(output.dir, "LOLA", "results_repName_LTR_hg19.rds"))


#Plot genomic regions in bubble plot
#subset result of interest
for(i in names(combined_analysis)){
dir.create(file.path(output.dir, i,"LOLA"), recursive=TRUE)

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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i,"LOLA", "results_repFamily_LTR_hg19_top20.pdf"), height=10)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA", "results_repFamily_LTR_hg19_topInf.pdf"), height=20)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA", "results_repFamily_LTR_hg19_Sign.pdf"), height=20)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA", "results_repClass_LTR_hg19_top20.pdf"), height=10)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA", "results_repClass_LTR_hg19_topInf.pdf"), height=20)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA", "results_repClass_LTR_hg19_Sign.pdf"), height=20)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA", "results_repName_LTR_hg19_top20.pdf"), height=10)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA","results_repName_LTR_hg19_top5.pdf"), height=10)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA","results_repName_LTR_hg19_topInf.pdf"), height=20)
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
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
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
pdf(file.path(output.dir, i, "LOLA","results_repName_LTR_hg19_Sign.pdf"), height=20)
print(g)
dev.off()

print(i)
}





#combined plot
anno_transcript <- anno[anno$type == "transcript",]
userUnisverse <- resize(anno_transcript,1)
#get degs
DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
UserSet <- lapply(DEG_results_list_sub, function(x){
#run enrichment analysis
resize(anno_transcript[anno_transcript$transcript_id %in%rownames(x[which(x$padj < alpha & x$log2FoldChange>lfc),]),],1)
results_repFamily_hg19= runLOLA(UserSet, userUnisverse, regionDB_repFamily_hg19, cores=4)
results_repClass_hg19= runLOLA(UserSet, userUnisverse, regionDB_repClass_hg19, cores=4)
results_repName_LTR_hg19= runLOLA(UserSet, userUnisverse, regionDB_repName_LTR_hg19, cores=4)

#prepare plotting
#results_repFamily_hg19
results <- results_repFamily_hg19
#idx_results_repFamily_hg19<- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repFamily_hg19<-head(unique(results$filename),1)
head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx_results_repFamily_hg19,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = rev(c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO")))
#prepare plotting
g_results_repFamily_hg19 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
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

#results_repClass_hg19
results <- results_repClass_hg19
#idx_results_repClass_hg19 <- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repClass_hg19 <-head(unique(results$filename),1)
results_sub <- results[results$filename %in% idx_results_repClass_hg19,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = rev(c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO")))
#prepare plotting
g_results_repClass_hg19 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
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


#results_repName_LTR_hg19
results <- results_repName_LTR_hg19
#idx_results_repName_LTR_hg19 <- unique(results[results$qValue < 0.05, ]$filename)
#idx_results_repName_LTR_hg19 <-head(unique(results$filename),5)
idx_results_repName_LTR_hg19 <-unique(results$filename[grep("LTR12",results$filename)])

results_sub <- results[results$filename %in% idx_results_repName_LTR_hg19,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = rev(c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO")))
#prepare plotting
g_results_repName_LTR_hg19 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+
scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
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

g <-cowplot::plot_grid(g_results_repFamily_hg19+coord_flip()+rremove("legend")+rremove("xlab"),
        g_results_repClass_hg19+coord_flip()+rremove("y.text")+rremove("legend")+rremove("xlab"),
        g_results_repName_LTR_hg19+coord_flip()+rremove("y.text")+rremove("xlab"),
        ncol=3,align = c("h"),#labels=c("Family", "Class", "Name"),
        rel_widths=c(length(idx_results_repFamily_hg19)+4, length(idx_results_repClass_hg19), length(idx_results_repName_LTR_hg19)))

ggsave(plot=g,file.path(output.dir,"LOLA", "repeats_combined_DEG_up.pdf"),height = 6, width = 10, useDingbats = FALSE)
