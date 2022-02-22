#plot coding potential score
#libraries
library(ggpubr)
library(data.table)
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(randomcoloR)
library(dendextend)
library(limma)
library(rafalib)
library(TxDb.Hsapiens.UCSC.hg19.knownGene )
library(seqinr)
#folder
base.dir<- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
pheno <- colData(vst)
dds <-readRDS(file.path(results.dir, "dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
fasta <- seqinr::read.fasta("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.fasta")
#housekeeper
hkt <- data.table::fread(file.path("/omics/groups/OE0219/internal/tinat/Encode/220128_deNovo_quantification_H1299ref_analysis/results/tables","MostStable.csv" ))

#get our transcript anno for housekeeper
anno_classi_known <- anno_classi[anno_classi$class_code =="=", ]
hkt_own <- anno_classi_known$transcript_id[sapply(strsplit(as.character(anno_classi_known$ref_id), ".", fixed=TRUE), "[",1) %in% hkt$"Ensembl ID" ]
length(hkt_own)

#get our transcript anno for non-coding
set.seed(42)
lincRNA <- unique(anno_original[which(anno_original$transcript_type=="lincRNA"),]$transcript_id)
lincRNA <- anno_classi_known[anno_classi_known$ref_id %in% lincRNA ,]$ref_id
lincRNA <-lincRNA[lincRNA %in%  anno_classi_known$qry_id]
lincRNA_sub <- sample(lincRNA, 100)

 #Take a look at design and annotation of samples
design(dds)
#Set specifications
cutoff <- 0.01 #set FDR cutoff
l2fc <- 2##set logfold2 cutoff

#look at overlap of upregulated degs and transcript class annotation
#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id

#combine deg with classcode information
DEG_results_list_anno <- lapply(DEG_results_list, function(x){
    x <- dplyr::left_join(x, anno_classi, by="transcript_id")
    x
})
DEG_results_list_anno_sub <- lapply(DEG_results_list_anno, function(x){
    x <- x[which(x$padj < cutoff & x$log2FoldChange >l2fc),]
    x <- split(x, x$class_code_simple)
    x
})
transcript_ofInterest <- lapply(DEG_results_list_anno_sub$DACandSB939_vs_DMSO, function(x){
    x$transcript_id
})
transcript_ofInterest$house_keeper <-hkt_own
transcript_ofInterest$lincRNA <- lincRNA_sub


#get housekeeping genes
#calculate protein coding potential for transcripts
#read fasta file into R
dir.create(file.path(PostDE.dir, "transcript_fasta"))
for(i in names(transcript_ofInterest)){
    write.fasta(fasta[transcript_ofInterest[[i]]], names=transcript_ofInterest[[i]],
    file.out=file.path(PostDE.dir, "transcript_fasta", paste0(i, ".fasta")) )
}



#read in data
results_cpc <- list.files(path=file.path(PostDE.dir, "transcript_fasta", "CPC2_output"), pattern="txt",full.names=TRUE)
results_cpc <- lapply(results_cpc, function(x){
    y <- fread(x)
    y$group <- sapply(strsplit(x, "/", fixed=TRUE),"[",13)
    y$group <- sapply(strsplit(y$group, ".", fixed=TRUE),"[",1)
    y$group <- sapply(strsplit(y$group, "cpc2_", fixed=TRUE),"[",2)
    y
    }
    )
results_cpc <- do.call("rbind", results_cpc)
 unique(results_cpc$group)
#plot coding propability
results_cpc$group <- gsub("^nonchimeric", "non-chimeric (novel)", results_cpc$group)
results_cpc$group <- gsub("^chimeric", "chimeric (novel)", results_cpc$group)

pal <- c("darkgray","gray", "orange", "tomato3")
names(pal)<- c("housekeeper","lincRNA",  "chimeric (novel)", "non-chimeric (novel)")
pdf(file.path(PostDE.dir, "transcript_fasta", "CPC2_output", "coding_probability.pdf"), height=5,width=5)
ggboxplot(results_cpc, x="group",y="coding_probability",ylab="Coding probability",
    fill="group", palette=pal, order=rev(c( "housekeeper", "lincRNA", "chimeric (novel)" , "non-chimeric (novel)" )))+
    rremove("legend") + 
    coord_flip() +
    rremove("ylab")
dev.off()
pdf(file.path(PostDE.dir, "transcript_fasta", "CPC2_output", "fickett_score.pdf"), height=5,width=5)
ggboxplot(results_cpc, x="group",y="Fickett_score",ylab="Fickett score",
    fill="group", palette=pal, order=rev(c( "housekeeper", "lincRNA", "chimeric (novel)" , "non-chimeric (novel)" )))+
    rremove("legend") + 
    coord_flip() +
    rremove("ylab")
dev.off()









#install rnamining 
#cd tools
#git clone https://gitlab.com/integrativebioinformatics/RNAmining.git
#cd RNAmining/
conda create -y --name RNAmining python>=3.8.0 pandas>=0.23.3 scikit-learn>=0.21.3 xgboost>=1.2.0 biopython>=1.78
conda activate RNAmining
#install dependencies
#conda install -c anacond pandas=0.23.3
#conda install -c anaconda scikit-learn=0.21.3
#conda install -c conda-forge xgboost=1.2.0
#conda install -c conda-forge biopython=1.78
#install rnamining
cd tools
git clone https://gitlab.com/integrativebioinformatics/RNAmining/
cd RNAmining
cd  volumes/rnamining-front/assets/scripts/

#run prediction on test data
head -n 100 /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.fasta > \
/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_test.fasta
mkdir /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/coding_prediction_test
python3 rnamining.py -f /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_test.fasta \
    -organism_name Homo_sapiens -prediction_type coding_prediction -output_folder /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/coding_prediction_test/

mkdir /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/coding_prediction
python3 rnamining.py -f /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.fasta \
    -organism_name Homo_sapiens -prediction_type coding_prediction -output_folder /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/coding_prediction/
