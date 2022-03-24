library(msa)
library(ggmsa)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome)
library(seqinr)
library(ape)
#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
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

#directories
base.dir <- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
dir.create(file.path(base.dir, "msa"))
#load data
repeats<- readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))
genome <- BSgenome.Hsapiens.UCSC.hg19 

#get repeat sequences
repeat_ltr12 <- repeats[grep("LTR12", repeats$repName),]
ltr12_sequences <- getSeq(genome,repeat_ltr12)
msa_default <- readRDS( file.path(base.dir, "msa", "msa_default.rds"))

#convert to align format
msa_default_align <- msaConvert(msa_default, type="seqinr::alignment")

#compute a distance matrix u
d <- dist.alignment(msa_default_align, "identity")
saveRDS(d, file.path(base.dir, "msa", "dist_alignment_identity.rds"))
ltr_Tree <- njs(d)
saveRDS(ltr_Tree, file.path(base.dir, "msa", "dist_alignment_identity_tree.rds"))
ltr_Tree <- readRDS(file.path(base.dir, "msa", "dist_alignment_identity_tree.rds"))

#plot tree
ltr_col <- brewer.pal(length(unique(as.character(repeat_ltr12[ltr_Tree$tip.label,]$repName))), "Accent")
names(ltr_col)<-unique(as.character(repeat_ltr12[ltr_Tree$tip.label,]$repName))
color_anno <- ltr_col[as.character(repeat_ltr12[ltr_Tree$tip.label,]$repName)]

pdf(file.path(base.dir, "msa", "dist_alignment_identity_tree.pdf"))
plot(ltr_Tree, main="Phylogenetic Tree of LTR12 Sequences", tip.color =color_anno)
dev.off()

