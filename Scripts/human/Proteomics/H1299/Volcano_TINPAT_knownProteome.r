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
library(UpSetR)

set.seed(42)

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
datasets.dir <- "/home/heyj/c010-datasets/Internal/COPD/enrichment_databases/"
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"

#Read in Data
#vst <- readRDS(file.path(results.dir, "vst.rds"))
#pheno <- colData(vst)
#dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
#tpm_mean_h1299 <- readRDS(file.path(results.dir, "tpm_meanTreatment_from_counts.rds"))

#load peptide data
#peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))

#load proteomics data
proteomics <- as.data.frame(data.table::fread(file.path("/omics/groups/OE0219/internal/tinat/proteomics/data/Limma_analysis_SBDAC_nascent_proteome.txt")))

#load aim and cta
ext.dir <- "/omics/groups/OE0219/internal/tinat/external_data/"
AIM <- as.data.frame(data.table::fread(file.path(ext.dir, "Li_etal", "ubiquitous_AIMS.txt")))
CTA <-  as.data.frame(data.table::fread(file.path(ext.dir,  "Cancer_Testis_Antigens.txt"), header=FALSE))

#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- factor(c("known",  "chimeric (novel)", "non-chimeric (novel)"), levels=c( "non-chimeric (novel)", "known",  "chimeric (novel)"))
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")


#subset peptides that originatee from our orf list
#peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
#select 2/3 candidates
#length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]$sequences))
#peptides_new_ORF_oi <- peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]

write.table(proteomics, file.path("/omics/groups/OE0219/internal/tinat/proteomics","proteomics_anno.csv"), sep="; ", col.names=TRUE, row.names=FALSE)
#define padj cutoff for plotting --> only up
alpha <- 0.05 #set FDR cutoff
lfc <- 0##set logfold2 cutoff
cutoff <- alpha
l2fc <- lfc

#Plot results in Volcano and MA plot {.tabset}
proteomics$padj <- (-(log10(proteomics$adj.P.Val )))
proteomics$col <- "not_significant"
proteomics$col <- ifelse(test = proteomics$adj.P.Val < alpha & abs(proteomics$logFC)>lfc, yes = "significant", 
                  no ="not_significant")
proteomics$SYMBOL  <- NA
proteomics[proteomics$Gene_names %in%CTA$V1, ]$SYMBOL <- proteomics[proteomics$Gene_names %in%CTA$V1, ]$Gene_names

  Volcanos <- ggplot(proteomics, aes(x=logFC, y=padj))+
    theme_pubr()+
    geom_point(data = subset(proteomics,col=="not_significant"), 
        fill = class_col[ "known"],color = class_col[ "known"], alpha = 0.2, size=2) +
    geom_point(data = subset(proteomics, col=="significant"), 
        fill = class_col["known"], color = class_col["known"],alpha = 1, size=2) +
    #ggrepel::geom_text_repel(aes(label = Label),
    #              box.padding   = 0.35, 
    #              point.padding = 0.5,
    #              segment.color = 'grey50') +
    geom_vline(xintercept=c(-lfc,lfc), linetype = 2) +
    geom_hline(yintercept=-log10(alpha), linetype = 2)+
    xlab("Protein expression change\n(log2 fold change)") + ylab("- log10(adj. P value)") +
    ggrepel::geom_text_repel(data = subset(proteomics, col=="significant"& !is.na(proteomics$SYMBOL) ),
    aes(label = SYMBOL),
    size = 3,
    seed = 42, 
    box.padding = 0.5,
    #max.overlaps = Inf,
    #arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .05,
    nudge_y = .1,
    color = "black"
  )
  
  ggsave(plot=Volcanos,file.path("/omics/groups/OE0219/internal/tinat/proteomics","volcano_knownProteome_Class_code_simple.pdf"),height = 5, width = 5, useDingbats = FALSE)
  ggsave(plot=Volcanos,file.path("/omics/groups/OE0219/internal/tinat/proteomics","volcano_knownProteome_Class_code_simple.png"),height = 5, width = 5, device="png")
