#comparison of peptides identified in patients and h1299
#libraries
library(ggpubr)
library(Biostrings)

#load data
peptides_new_AML_MHCI <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/integration/peptidomics_AML/comparison_gene_expression/MHCI", "peptides_list_new.rds"))
peptides_new_H1299 <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression", "peptides_list_new.rds"))
lift <- readRDS("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly//H1299_comparison_annotatedVersion/liftAnno_celllinePanel_h1299.rds")


#subset h1299 peptides
peptides_new_H1299_ORF <- peptides_new_H1299[peptides_new_H1299$Species =="ORFs",]
peptides_new_H1299_ORF_oi <- peptides_new_H1299_ORF[peptides_new_H1299_ORF$DAC_SB>60 & peptides_new_H1299_ORF$DMSO ==0,]
nrow(peptides_new_H1299_ORF_oi)
length(unique(peptides_new_H1299_ORF_oi$sequence))
#subset aml peptides
peptides_new_AML_MHCI_ORF <- peptides_new_AML_MHCI[peptides_new_AML_MHCI$Species =="ORFs",]
treat_cand_seq <- peptides_new_AML_MHCI_ORF[peptides_new_AML_MHCI_ORF$"Timepoint [h]">47 ,]$Sequence
untreat_cand_seq <- peptides_new_AML_MHCI_ORF[peptides_new_AML_MHCI_ORF$"Timepoint [h]"==0 ,]$Sequence
sel_cand <- treat_cand_seq[!treat_cand_seq %in% untreat_cand_seq]
peptides_new_AML_MHCI_ORF_oi <- peptides_new_AML_MHCI_ORF[peptides_new_AML_MHCI_ORF$Sequence %in% sel_cand,]
nrow(peptides_new_AML_MHCI_ORF_oi)
length(unique(peptides_new_AML_MHCI_ORF_oi$Sequence))

#combine peptide sequences with lift information
lift_equal <- lift[lift$class_code=="=",]
peptides_new_H1299_ORF_oi$H1299_anno <- peptides_new_H1299_ORF_oi$transcript_id
peptides_new_AML_MHCI_ORF_oi$Cellline_anno <- peptides_new_AML_MHCI_ORF_oi$transcript_id
peptides_new_AML_MHCI_ORF_oi <- dplyr::left_join(peptides_new_AML_MHCI_ORF_oi, lift_equal, by="Cellline_anno")

#check stats
#aml perspective
table(peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,useNA = c("always"))
peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$H1299_anno
unique(peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$H1299_anno)
peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$Cellline_anno
unique(peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$Cellline_anno)
peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$accession_new
unique(peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$accession_new)
unique(peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$Donor)
peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$Sequence
unique(peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$Sequence)
unique(peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$accession)
write.table(peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,], 
    file.path("/omics/groups/OE0219/internal/tinat/integration/", "H1299_AML_peptide_comparison", "AML_peptides_MHCI_overlapH1299_transcript.txt"), sep="\t",quote=FALSE, row.names=FALSE)


#h1299perspective
table(peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,useNA = c("always"))
peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,]$H1299_anno
unique(peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,]$H1299_anno)
peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,]$accession_new
unique(peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,]$accession)
peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,]$sequence
unique(peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,]$sequence)
write.table(peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,], 
    file.path("/omics/groups/OE0219/internal/tinat/integration/", "H1299_AML_peptide_comparison", "H1299_peptides_overlapAMLMHCI_transcript.txt"), sep="\t",quote=FALSE, row.names=FALSE)

#plot overlap
transcript_list  <- list(AML=peptides_new_AML_MHCI_ORF_oi$H1299_anno, H1299=peptides_new_H1299_ORF_oi$H1299_anno)
dir.create(file.path("/omics/groups/OE0219/internal/tinat/integration/", "H1299_AML_peptide_comparison"))
pdf(file.path(file.path("/omics/groups/OE0219/internal/tinat/integration/", "H1299_AML_peptide_comparison"), 
    "Venn_transcriptPeptides.pdf"))
ggVennDiagram::ggVennDiagram(transcript_list, label_alpha = 0, category.names=names(transcript_list))+
    rremove("legend") + 
    theme(text = element_text(size = 8)) #+ 
    #scale_fill_gradient(low="white",high = "#EF8A62")
dev.off() 

#get orf sequences
pepteide_sub_induced_AML <- readAAStringSet("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACSBvsDMSO_forJens_novel.fa")
pepteide_sub_induced_H1299 <- readAAStringSet( "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa")

#subset orfs
#for AML
temp <- sapply(strsplit(names(pepteide_sub_induced_AML), " type",fixed=TRUE), function(x)x[1]) 
pepteide_sub_induced_AML_sub <- pepteide_sub_induced_AML[temp %in% peptides_new_AML_MHCI_ORF_oi[peptides_new_AML_MHCI_ORF_oi$H1299_anno %in% peptides_new_H1299_ORF_oi$H1299_anno,]$accession_new,]
#for H1299

temp <- sapply(strsplit(names(pepteide_sub_induced_H1299), " type",fixed=TRUE), function(x)x[1]) 
pepteide_sub_induced_H1299_sub <- pepteide_sub_induced_H1299[temp %in% peptides_new_H1299_ORF_oi[peptides_new_H1299_ORF_oi$H1299_anno %in% peptides_new_AML_MHCI_ORF_oi$H1299_anno,]$accession_new,]

#look at overlap
length(pepteide_sub_induced_AML_sub)
length(pepteide_sub_induced_H1299_sub)
table(pepteide_sub_induced_AML_sub %in% pepteide_sub_induced_H1299_sub)
pepteide_sub_induced_AML_sub[pepteide_sub_induced_AML_sub %in% pepteide_sub_induced_H1299_sub]
writeXStringSet(pepteide_sub_induced_AML_sub[pepteide_sub_induced_AML_sub %in% pepteide_sub_induced_H1299_sub],  
    file.path("/omics/groups/OE0219/internal/tinat/integration/", "H1299_AML_peptide_comparison", "commonORFs_AML.fa"))
table(pepteide_sub_induced_H1299_sub %in% pepteide_sub_induced_AML_sub)
pepteide_sub_induced_H1299_sub[pepteide_sub_induced_H1299_sub %in%pepteide_sub_induced_AML_sub]
writeXStringSet(pepteide_sub_induced_H1299_sub[pepteide_sub_induced_H1299_sub %in%pepteide_sub_induced_AML_sub],  
    file.path("/omics/groups/OE0219/internal/tinat/integration/", "H1299_AML_peptide_comparison", "commonORFs_H1299.fa"))


#compare on orf level directly
#subset orfs
#for AML
temp <- sapply(strsplit(names(pepteide_sub_induced_AML), " type",fixed=TRUE), function(x)x[1]) 
length(unique(peptides_new_AML_MHCI_ORF_oi$accession))
length(unique(peptides_new_AML_MHCI_ORF_oi$Sequence))
pepteide_sub_induced_AML_sub <- pepteide_sub_induced_AML[temp %in% peptides_new_AML_MHCI_ORF_oi$accession_new,]
length(pepteide_sub_induced_AML_sub)
#subset
pepteide_sub_induced_AML_sub_unique <- unique(pepteide_sub_induced_AML_sub)
length(pepteide_sub_induced_AML_sub_unique)
#for H1299
temp <- sapply(strsplit(names(pepteide_sub_induced_H1299), " type",fixed=TRUE), function(x)x[1]) 
length(unique(peptides_new_H1299_ORF_oi$accession))
length(unique(peptides_new_H1299_ORF_oi$sequence))
pepteide_sub_induced_H1299_sub <- pepteide_sub_induced_H1299[temp %in% peptides_new_H1299_ORF_oi$accession_new,]
length(pepteide_sub_induced_H1299_sub)
#subset
pepteide_sub_induced_H1299_sub_unique <- unique(pepteide_sub_induced_H1299_sub)
length(pepteide_sub_induced_H1299_sub_unique)

pepteide_sub_induced_AML_sub_unique[pepteide_sub_induced_AML_sub_unique %in% pepteide_sub_induced_H1299_sub_unique]
pepteide_sub_induced_H1299_sub_unique[pepteide_sub_induced_H1299_sub_unique %in% pepteide_sub_induced_AML_sub]


