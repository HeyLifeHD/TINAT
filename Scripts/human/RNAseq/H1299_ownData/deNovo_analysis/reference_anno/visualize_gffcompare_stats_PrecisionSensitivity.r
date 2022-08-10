#visualize gffcompare stats

#library
library(ggpubr)

#directories
output.dir <- "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/visualizations_custom"
dir.create(output.dir)

#load input
stats <-  as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.stats"))

#extract sensitivity precision matrix
stats1 <- stats[1:6,]
level<-sapply(strsplit(stats1[,1],":"),"[",1)
sensitivity<- sapply(strsplit(stats1[,1],":"),"[",2)
sensitivity <- gsub(" ", "",sensitivity)
sensitivity <- as.numeric(gsub("\"", "",sensitivity))
precision<- as.numeric(stats1[,2])
stats1 <-data.frame(level=as.character(level), sensitivity=sensitivity, precision=precision)

#prepare for plotting
plot <- reshape2::melt(stats1)
pdf(file.path(output.dir, "precision_sensitivity_barplot.pdf"))
ggbarplot(plot,x="level", y="value", fill="variable",  label = TRUE, 
    ylab="Fraction of [%]", ylim=c(0,100), order=rev(stats1$level),
    position = position_dodge(0.9), palette=c("#E7B800", "#FC4E07")) + 
    rotate() +
    rremove("legend.title") +
    rremove("ylab")
dev.off()

#plot novel missed and total
stats1 <- stats[grep("Novel", stats[,1]),]
level  <-sapply(strsplit(stats1[,1],":"),"[",1)
level <- gsub("Novel ","",level)
level <- stringr::str_to_title(level)
novel<- sapply(strsplit(stats1[,1],":"),"[",2)
novel<- sapply(strsplit(novel,"/"),"[",1)
novel <- gsub(" ","", novel)
total<- sapply(strsplit(stats1[,1],":"),"[",2)
total<- sapply(strsplit(total,"/", fixed=TRUE),"[",2)
total<- sapply(strsplit(total,"\t", fixed=TRUE),"[",1)
stats1 <- stats[grep("Missed", stats[,1]),]
missed<- sapply(strsplit(stats1[,1],":"),"[",2)
missed<- sapply(strsplit(missed,"/"),"[",1)
missed <- gsub(" ","", missed)
stats1 <-data.frame(level=as.character(level), novel=as.integer(novel), 
    missed=as.integer(missed), total=as.integer(total))


#prepare for plotting
plot <- reshape2::melt(stats1, id.vars= c("level"))
pdf(file.path(output.dir, "novel_missed_loci_barplot.pdf"))
ggbarplot(plot,x="level", y="value", fill="variable",  label = TRUE, 
    ylab="Number of", ylim=c(0,100), order=rev(stats1$level),
    position = position_dodge(0.9), palette=c("#00AFBB", "#E7B800", "#FC4E07")) + 
    rotate() +
    rremove("legend.title") +
    rremove("ylab")
dev.off()
