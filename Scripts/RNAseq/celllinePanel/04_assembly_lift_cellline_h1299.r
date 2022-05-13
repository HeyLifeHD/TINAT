
#gffcompare 
mkdir /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison 
conda activate gffcompare
cp  /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison/reference.gtf
cp /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/mergedTranscripts.gtf /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison/input.gtf
gffcompare -r  /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison/reference.gtf \
-o /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison/gffCompare /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison/input.gtf
#sort gff compare
/home/heyj/tools/gff3sort/gff3sort.pl /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison/gffCompare.annotated.gtf  > /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison/gffCompare.annotated.sorted.gtf

#gffcompare of annotated version
mkdir /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion
cp  /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion/reference.gtf
cp /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.gtf /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion/input.gtf
conda activate gffcompare
gffcompare -r /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion/reference.gtf \
-o /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion/gffCompare /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion/input.gtf
#sort gff compare
/home/heyj/tools/gff3sort/gff3sort.pl /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion/gffCompare.annotated.gtf  > /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/H1299_comparison_annotatedVersion/gffCompare.annotated.sorted.gtf

#generate dictionary for lift over of cell line to h1299 anno
#comparison
base.dir<- "/omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification_analysis/"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"

#Read in Data
vst <- readRDS(file.path(results.dir, "vst.rds"))
anno <- readRDS("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.mergedTranscripts.gtf.tmap"))

anno_new <- rtracklayer::import.gff2("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly//H1299_comparison_annotatedVersion/gffCompare.annotated.sorted.gtf")
anno_h1299 <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")

#compare rownames
transcripts <- rownames(vst)
length(transcripts)
anno_new_transcript <- anno_new[anno_new$type=="transcript",]
length(anno_new_transcript)
table(transcripts %in% anno_new_transcript$ref_gene_id)
table(transcripts %in% anno_new_transcript$cmp_ref)
table(transcripts %in% anno_new_transcript$transcript_id)
#generate dictionary
lift <- data.frame(Cellline_anno=as.character(anno_new_transcript$transcript_id), H1299_anno=as.character(anno_new_transcript$cmp_ref), class_code=as.character(anno_new_transcript$class_code))
lift$Cellline_anno <- as.character(lift$Cellline_anno)
lift$H1299_anno <- as.character(lift$H1299_anno)
saveRDS(lift, "/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly//H1299_comparison_annotatedVersion/liftAnno_celllinePanel_h1299.rds")
lift <- readRDS("/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly//H1299_comparison_annotatedVersion/liftAnno_celllinePanel_h1299.rds")
write.table(lift,"/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly//H1299_comparison_annotatedVersion/liftAnno_celllinePanel_h1299.tsv",
    sep="\t", row.names=FALSE, col.names=TRUE )

#see how many peptide candidate transcripts have a perfect match
#load peptide data
peptides_new <- readRDS(file.path(output.dir, "peptides_list_new.rds"))
#subset peptides that originatee from our orf list
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
#select 2/3 candidates
length(unique(peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]$sequences))
peptides_new_ORF_oi <- peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]


#subset equal match transcripts
lift_exact <- lift[lift$class_cod=="=",]
dim(lift_exact)
nrow(peptides_new_ORF_oi[which(peptides_new_ORF_oi$transcript_id %in% lift_exact$H1299_anno),])
nrow(peptides_new_ORF_oi)

length(unique(peptides_new_ORF_oi[which(peptides_new_ORF_oi$transcript_id %in% lift_exact$H1299_anno),]$sequences))
length(unique(peptides_new_ORF_oi$sequences))

unique(peptides_new_ORF_oi[!peptides_new_ORF_oi$transcript_id %in% lift_exact$H1299_anno,]$transcript_id)

# table(transcripts %in% lift$Cellline_anno)
# str(lift)

#  lift[which(lift$H1299_anno == "MSTRG.23152.1" & lift$class_code=="="),]
#  anno_h1299[anno_h1299$transcript_id == "MSTRG.23152.1" ,]
# saveRDS()


# problematic:: MSTRG.23152.1; MSTRG.29191.1; MSTRG.29191.2

# anno_new_transcript[anno_new_transcript %over% anno_h1299[anno_h1299$transcript_id == "MSTRG.23152.1" ,], ]
# anno_h1299 <-  readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
# anno_h1299_transcript <- anno_h1299[anno_h1299$type=="transcript",]
# anno_new_transcript[anno_new_transcript %over% anno_h1299_transcript[anno_h1299_transcript$transcript_id =="MSTRG.23152.1" ,],]
#  chr3 114530457-114531041 
# ENST00000357258.8_3