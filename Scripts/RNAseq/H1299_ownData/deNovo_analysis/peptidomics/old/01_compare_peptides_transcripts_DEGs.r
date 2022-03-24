#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)
#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"

output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"
dir.create(output.dir)
#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))

#peptidomics list
peptides <- as.data.frame(readxl::read_excel("/omics/groups/OE0219/internal/tinat/integration/peptidomics/data/220301_Peptide_Results_new.xlsx", sheet=1))

#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id

#prepare list to have a row for each accession
peptides_new <- list()
for(i in 1:length(strsplit(peptides$Accession,";", fixed=TRUE))){
  if(length(strsplit(peptides$Accession,";", fixed=TRUE)[[i]])==1){
    peptides_new[[i]]<- peptides[i,]
    peptides_new[[i]]$accession_new <-  peptides[i,]$Accession
    peptides_new[[i]]$origin<-  1
    peptides_new[[i]]$total_origins <- 1
  } else if(length(strsplit(peptides$Accession,";", fixed=TRUE)[[i]])>1){
    peptides_new[[i]]<- peptides[rep(i,length(strsplit(peptides$Accession,";", fixed=TRUE)[[i]])),]                    
    peptides_new[[i]]$accession_new<-  strsplit(peptides$Accession,";", fixed=TRUE)[[i]]
    peptides_new[[i]]$origin<-  1:length(strsplit(peptides$Accession,";", fixed=TRUE)[[i]])
    peptides_new[[i]]$total_origins <- length(strsplit(peptides$Accession,";", fixed=TRUE)[[i]])
  }
  #if(length(strsplit(peptides_new[[i]]$accession_new,"ORF_", fixed=TRUE))>1){
  peptides_new[[i]]$transcript_id <- sapply(strsplit(peptides_new[[i]]$accession_new,"ORF_", fixed=TRUE),"[",2)
  peptides_new[[i]]$transcript_id <- sapply(strsplit( peptides_new[[i]]$transcript_id,".p", fixed=TRUE),"[",1)
  #} 
}
peptides_new <- do.call("rbind", peptides_new)
#get gene id for protein names
#library(AnnotationDbi)
#library(org.Hs.eg.db)
#mapIds(org.Hs.eg.db, keys= peptides_new$accession_new , keytype ="UNIPROT", column = "SYMBOL", multiVals = "first" )

#combine annotation with peptides
peptides_new <- dplyr::left_join(peptides_new, anno_classi, by="transcript_id")
saveRDS(peptides_new, file.path(output.dir, "peptides_list_new.rds"))
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))
write.table(peptides_new, file.path(output.dir, "peptides_list_new.tsv"), sep="\t",quote=FALSE, row.names=FALSE)


nrow(peptides[peptides$DAC_SB==100 & peptides$DMSO ==0,])
length(unique(peptides[peptides$DAC_SB==100 & peptides$DMSO ==0 ,]$sequences))
peptides_ORF <- peptides[peptides$Species =="ORFs",]

length(unique(peptides_ORF[peptides_ORF$DAC_SB==100 & peptides_ORF$DMSO ==0,]$sequences))

#subset peptides that originatee from our orf list
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB==100 & peptides_new_ORF$DMSO ==0,]$sequences))
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB==100 & peptides_new_ORF$DMSO ==0,]$transcript_id))
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]$transcript_id))
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]$transcript_id))

#look at deg comparison
#define padj cutoff for plotting --> only up
cutoff <- alpha
l2fc <- lfc
DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
genes2plot <- lapply(DEG_results_list_sub, function(x){
  x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
  x <- rownames(x)
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]
lapply(genes2plot, function(x)length(x))

#get peptides unique to DAC_SB
peptides_new_sub <- peptides_new[peptides_new$DAC_SB>0 & peptides_new$DMSO ==0,]
compar <- c(genes2plot,list(uniqueDACSB_peptides = unique(peptides_new_sub$transcript_id)))
#plot
#col <- c("#00AFBB", "#E7B800", "#FC4E07", "orange", "tomato3")#"darkgray"
#names(col)<- c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO", "chimeric (novel)", "non-chimeric (novel)")# "known", 
#col <- c( "#FC4E07","#00AFBB", "tomato3", "#FC4E07", "orange","darkgray")
pdf(file.path(output.dir, "Upset_upDEG_DACSB_unique_peptides.pdf"), width=10, height=5)
UpSetR::upset(UpSetR::fromList(compar),  order.by = "freq", set_size.show = TRUE ,nsets=6)
dev.off()
#same but in all three dac sb replicates
peptides_new_sub <- peptides_new_ORF[peptides_new$DAC_SB==100 & peptides_new$DMSO ==0,]
compar <- c(genes2plot,list(uniqueDACSB_peptides = unique(peptides_new_sub$transcript_id)))
#plot
pdf(file.path(output.dir, "Upset_upDEG_DACSB_unique_peptides_DACSB100.pdf"), width=10, height=5)
UpSetR::upset(UpSetR::fromList(compar),  order.by = "freq", set_size.show = TRUE ,nsets=6,)
dev.off()


#same for peptides from our orf list
#get peptides unique to DAC_SB
peptides_new_sub <- peptides_new_ORF[peptides_new_ORF$DAC_SB>0 & peptides_new_ORF$DMSO ==0,]
compar <- c(genes2plot,list(uniqueDACSB_peptides = unique(peptides_new_sub$transcript_id)))
#plot
pdf(file.path(output.dir, "Upset_upDEG_DACSB_unique_peptides_ourORF.pdf"), width=10, height=5)
UpSetR::upset(UpSetR::fromList(compar),  order.by = "freq", set_size.show = TRUE ,nsets=6)
dev.off()
###to do --> modify source code to not scale data and only modify axis

#same but in all three dac sb replicates
peptides_new_sub <- peptides_new_ORF[peptides_new_ORF$DAC_SB==100 & peptides_new_ORF$DMSO ==0,]
compar <- c(genes2plot,list(uniqueDACSB_peptides = unique(peptides_new_sub$transcript_id)))
#plot
pdf(file.path(output.dir, "Upset_upDEG_DACSB_unique_peptides_ourORF_DACSB100.pdf"), width=10, height=5)
UpSetR::upset(UpSetR::fromList(compar),  order.by = "freq", set_size.show = TRUE ,nsets=6#,intersections=
  #list(list("uniqueDACSB_peptides", "DACandSB939_vs_DMSO"),list("uniqueDACSB_peptides","SB939_vs_DMSO"),list("uniqueDACSB_peptides","DAC_vs_DMSO" ))
  )
dev.off()
pdf(file.path(output.dir, "Upset_upDEG_DACSB_unique_peptides_ourORF_DACSB100_selectedIntersection.pdf"), width=10, height=5)
UpSetR::upset(UpSetR::fromList(compar),  order.by = "freq", set_size.show = TRUE ,nsets=6,intersections=
  list(list("uniqueDACSB_peptides", "DAC_vs_DMSO","DACandSB939_vs_DMSO"),list("uniqueDACSB_peptides","DACandSB939_vs_DMSO"),list("uniqueDACSB_peptides","DACandSB939_vs_DMSO","SB939_vs_DMSO"),list("uniqueDACSB_peptides","DAC_vs_DMSO","DACandSB939_vs_DMSO","SB939_vs_DMSO" ))
  )
dev.off()

#same for more conditions
compar<- c(genes2plot,list(DACSB_peptides_100 = unique(peptides_new_ORF[peptides_new_ORF$DAC_SB==100 & peptides_new_ORF$DMSO ==0,]$transcript_id)),
  list(DACSB_peptides_66 = unique( peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DAC_SB<70  & peptides_new_ORF$DMSO ==0,]$transcript_id)),
    list(DACSB_peptides_33 = unique( peptides_new_ORF[peptides_new_ORF$DAC_SB>30 & peptides_new_ORF$DAC_SB<40  & peptides_new_ORF$DMSO ==0,]$transcript_id)))
 pdf(file.path(output.dir, "Upset_upDEG_DACSB_unique_peptides_ourORF_DACSB.pdf"), width=10, height=5)
UpSetR::upset(UpSetR::fromList(compar),  order.by = "freq", set_size.show = TRUE ,nsets=6#,intersections=
  #list(list("uniqueDACSB_peptides", "DACandSB939_vs_DMSO"),list("uniqueDACSB_peptides","SB939_vs_DMSO"),list("uniqueDACSB_peptides","DAC_vs_DMSO" ))
  )
dev.off()

#merge with deg list and plot volcano
library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
DEG_results_list_plot <- DEG_results_list
DEG_results_list_plot <- lapply(DEG_results_list_plot, function(x){
  x <- dplyr::left_join(x, peptides_new_sub, by="transcript_id")
  x$col <- "gray"
  x$col <- ifelse(test = x$padj < alpha & x$log2FoldChange>lfc, yes =  col[1], 
                  no = ifelse(test = x$padj < alpha  & x$log2FoldChange<(-lfc), yes =  col[3], "gray"))
  x$col <- ifelse(is.na(x$sequences), x$col, "black")
  temp <- c(head(x[x$col== col[1],]$symbol, 5), head(x[x$col== col[3],]$symbol, 5))
  x$SYMBOL <- NA
  x[x$symbol %in% temp,]$SYMBOL <- x[x$symbol %in% temp,]$symbol
  x$pvalue <-(-(log10(x$pvalue)))
  x$padj <- (-(log10(x$padj)))
  x
})

#Volcano Plot
Volcanos<- list() 
for (i in names(DEG_results_list_plot)){
  dir.create(file.path(output.dir,i))
  Volcanos[[i]] <- ggplot(DEG_results_list_plot[[i]], aes(x=log2FoldChange, y=padj))+
    theme_pubr()+
    geom_point(data = subset(DEG_results_list_plot[[i]], col!="black"), 
             color = subset(DEG_results_list_plot[[i]], col!="black")$col, alpha = 0.5, size=.6) +
    geom_point(data = subset(DEG_results_list_plot[[i]], col=="black"), 
             color = "black", alpha =1, size=.6) + 
    geom_vline(xintercept=c(-lfc,lfc), linetype = 2) +
    geom_hline(yintercept=-log10(alpha), linetype = 2)+
     xlab("Gene expression change\n(log2 fold change)") + ylab("- log10(P.adjusted)")
  
  ggsave(plot=Volcanos[[i]],file.path(output.dir,i,"volcano_highlight_DACSB_uniquePeptides.pdf"),height = 4, width = 4, useDingbats = FALSE)
  ggsave(plot=Volcanos[[i]],file.path(output.dir,i,"volcano_highlight_DACSB_uniquePeptides.png"),height = 4, width = 4, device="png")
  print(i)
}

#MA plot
MAs<- list() 
for (i in names(DEG_results_list_plot)){
  dir.create(file.path(output.dir,i))
  MAs[[i]] <- ggplot(DEG_results_list_plot[[i]], aes(x=baseMean, y=log2FoldChange))+
    theme_pubr()+
      geom_point(data = subset(DEG_results_list_plot[[i]], col!="black"), 
             color = subset(DEG_results_list_plot[[i]], col!="black")$col, alpha = 0.5, size=.6) +
      geom_point(data = subset(DEG_results_list_plot[[i]], col=="black"), 
             color = "black", alpha =1, size=.6) + 
      geom_hline(yintercept=c(-lfc,lfc), linetype = 2)+
      xlab("Average mean expression [log10]") + ylab("Gene expression change\n(log2 fold change))")+
      scale_x_continuous(trans='log10')
  ggsave(plot=MAs[[i]],file.path(output.dir,i,"MA_highlight_DACSB_uniquePeptides.pdf"),height = 4, width = 4, useDingbats = FALSE)
  ggsave(plot=MAs[[i]],file.path(output.dir,i,"MA_highlight_DACSB_uniquePeptides.png"),height = 4, width = 4,device="png")
  print(i)
}


#plot class code simple of transcripts of upregualted degs with erg distance 
#get transcript and repeat annotation
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
anno_classi$transcript_id <- anno_classi$qry_id 

DEG_results_list_anno_sub_DACSBsep<- lapply(DEG_results_list, function(x){
  x <- dplyr::left_join(x, anno_classi, by="transcript_id")
  x$ERE_TSS_anno <- NA
  x$ERE_TSS_anno <- ifelse(x$dist_nearest_repeat == 0, x$nearest_repeat_repClass, "no ERE")
  x$ERE_TSS_anno <- ifelse(x$ERE_TSS_anno %in% c("LTR", "LINE", "SINE", "no ERE"), x$ERE_TSS_anno , "other")
  y <- list(DAC_SB100 = x[x$transcript_id %in% unique(peptides_new_ORF[peptides_new_ORF$DAC_SB==100 & peptides_new_ORF$DMSO ==0,]$transcript_id),],
    DAC_SB66 = x[x$transcript_id %in%  unique( peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DAC_SB<70  & peptides_new_ORF$DMSO ==0,]$transcript_id),],
     DAC_SB33 = x[x$transcript_id %in% unique( peptides_new_ORF[peptides_new_ORF$DAC_SB>30 & peptides_new_ORF$DAC_SB<40  & peptides_new_ORF$DMSO ==0,]$transcript_id),])
  y
})

#get class code statistics
class_code_simple_stat <- lapply(DEG_results_list_anno_sub_DACSBsep,function(x){
  x <- lapply(x, function(y){
    y <- as.data.frame(table(y$class_code_simple))
    y})
  x
})

#plot
pal <- c("darkgray", "orange", "tomato3")
names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
pie <- list()
for (i in names(class_code_simple_stat)){
    pie[[i]]<- list()
    for(j in names(class_code_simple_stat[[i]])){
        labs <- paste0(class_code_simple_stat[[i]][[j]]$Freq)
        pie[[i]][[j]] <- ggpie(class_code_simple_stat[[i]][[j]], "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=pal,
                fill = "Var1", color = "black",title=j)+rremove("legend.title")
    }
}


pdf(file.path(output.dir, "Pie_DAC_SB_peptide_origin_Class_code_simple_ORFs_differentOrder.pdf"), width=4, height=10)
    ggarrange(pie$"DACandSB939_vs_DMSO"$"DAC_SB33",  pie$"DACandSB939_vs_DMSO"$"DAC_SB66", pie$"DACandSB939_vs_DMSO"$"DAC_SB100",  
        common.legend = TRUE,
        labels =c("1/3 DAC + SB939 unique", "2/3 DAC + SB939 unique", "3/3 DAC + SB939 unique"),
        ncol = 1, nrow = 3
        )
dev.off()


#same for ere stats
ERE_simple_stat <- lapply(DEG_results_list_anno_sub_DACSBsep,function(x){
  lapply(x, function(y)as.data.frame(table(y$ERE_TSS_anno)))
})


#plot
pal <- c("deeppink3", "red4", "midnightblue", "gold1", "salmon4")
names(pal)<-  levels(ERE_simple_stat$DACandSB939_vs_DMSO$DAC_SB66$Var1)  

pie <- list()
for (i in names(ERE_simple_stat)){
    pie[[i]]<- list()
    for(j in names(ERE_simple_stat[[i]])){
        labs <- paste0(ERE_simple_stat[[i]][[j]]$Freq)
        pie[[i]][[j]] <- ggpie(ERE_simple_stat[[i]][[j]], "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=pal,
                fill = "Var1", color = "black",title=j)+rremove("legend.title")
    }
}


pdf(file.path(output.dir, "Pie_DAC_SB_peptide_origin_EREanno_simple_ORFs_differentOrder.pdf"), width=4, height=10)
    ggarrange(pie$"DACandSB939_vs_DMSO"$"DAC_SB33",  pie$"DACandSB939_vs_DMSO"$"DAC_SB66", pie$"DACandSB939_vs_DMSO"$"DAC_SB100",  
        common.legend = TRUE,
        labels = c("1/3 DAC + SB939 unique", "2/3 DAC + SB939 unique", "3/3 DAC + SB939 unique"),
        ncol = 1, nrow = 3
        )
dev.off()


#by DEG upregulated
#[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
DEG_results_list_anno_sub_DACSBsep_sub<- lapply(DEG_results_list_anno_sub_DACSBsep, function(x){
  lapply(x, function(y){
    y[which(y$padj < cutoff & y$log2FoldChange >l2fc),]
  })
})

#get class code statistics
class_code_simple_stat <- lapply(DEG_results_list_anno_sub_DACSBsep_sub,function(x){
  x <- lapply(x, function(y){
    y <- as.data.frame(table(y$class_code_simple))
    y})
  x
})

#plot
pal <- c("darkgray", "orange", "tomato3")
names(pal)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
pie <- list()
for (i in names(class_code_simple_stat)){
    pie[[i]]<- list()
    for(j in names(class_code_simple_stat[[i]])){
        labs <- paste0(class_code_simple_stat[[i]][[j]]$Freq)
        pie[[i]][[j]] <- ggpie(class_code_simple_stat[[i]][[j]], "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=pal,
                fill = "Var1", color = "black",title=j)+rremove("legend.title")
    }
}


pdf(file.path(output.dir, "Pie_DEGup_DAC_SB_peptide_origin_Class_code_simple_ORFs_differentOrder.pdf"), width=10, height=10)
    ggarrange(pie$"DAC_vs_DMSO"$"DAC_SB33",  pie$"DAC_vs_DMSO"$"DAC_SB66", pie$"DAC_vs_DMSO"$"DAC_SB100",  
    pie$"SB939_vs_DMSO"$"DAC_SB33",  pie$"SB939_vs_DMSO"$"DAC_SB66", pie$"SB939_vs_DMSO"$"DAC_SB100",  
    pie$"DACandSB939_vs_DMSO"$"DAC_SB33",  pie$"DACandSB939_vs_DMSO"$"DAC_SB66", pie$"DACandSB939_vs_DMSO"$"DAC_SB100",  
        common.legend = TRUE,
        labels = names(pie[[1]]),
        ncol = 3, nrow = 3
        )
dev.off()


#same for ere stats
ERE_simple_stat <- lapply(DEG_results_list_anno_sub_DACSBsep,function(x){
  lapply(x, function(y)as.data.frame(table(y$ERE_TSS_anno)))
})


#plot
pal <- c("deeppink3", "red4", "midnightblue", "gold1", "salmon4")
names(pal)<-  levels(ERE_simple_stat$DACandSB939_vs_DMSO$DAC_SB66$Var1)  

pie <- list()
for (i in names(ERE_simple_stat)){
    pie[[i]]<- list()
    for(j in names(ERE_simple_stat[[i]])){
        labs <- paste0(ERE_simple_stat[[i]][[j]]$Freq)
        pie[[i]][[j]] <- ggpie(ERE_simple_stat[[i]][[j]], "Freq", #label = labs,
                lab.pos = "out", lab.font = "white",label.size = 16, palette=pal,
                fill = "Var1", color = "black",title=j)+rremove("legend.title")
    }
}


pdf(file.path(output.dir, "Pie_DEGup_DAC_SB_peptide_origin_EREanno_simple_ORFs_differentOrder.pdf"), width=10, height=10)
    ggarrange(pie$"DAC_vs_DMSO"$"DAC_SB33",  pie$"DAC_vs_DMSO"$"DAC_SB66", pie$"DAC_vs_DMSO"$"DAC_SB100",  
    pie$"SB939_vs_DMSO"$"DAC_SB33",  pie$"SB939_vs_DMSO"$"DAC_SB66", pie$"SB939_vs_DMSO"$"DAC_SB100", 
    pie$"DACandSB939_vs_DMSO"$"DAC_SB33",  pie$"DACandSB939_vs_DMSO"$"DAC_SB66", pie$"DACandSB939_vs_DMSO"$"DAC_SB100", 
        common.legend = TRUE,
        labels = names(pie[[1]]),
        ncol = 3, nrow = 3
        )
dev.off()