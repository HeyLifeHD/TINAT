library(msa)
library(ggmsa)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome)
#directories
base.dir <- "/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/"
dir.create(file.path(base.dir, "msa"))
#load data
repeats<- readRDS(file.path("c010-datasets/Internal/LTR/guideDesign/","Workspaces_David","repeats.RDS"))
genome <- BSgenome.Hsapiens.UCSC.hg19 

#get repeat sequences
repeat_ltr12 <- repeats[grep("LTR12", repeats$repName),]
ltr12_sequences <- getSeq(genome,repeat_ltr12)

#run msa with default sequences
msa_default <- msa(ltr12_sequences, method="ClustalW", verbose=TRUE)
saveRDS(msa_default, file.path(base.dir, "msa", "msa_default.rds"))
msa_default <- readRDS( file.path(base.dir, "msa", "msa_default.rds"))
writeXStringSet(msa_default@unmasked,file.path(base.dir, "msa", "msa_default_unmasked.fasta"))

#get consensus matrices
conMat <- consensusMatrix(msa_default@unmasked)

#plot pretty msa
temp <-ggmsa(msa_default@unmasked, seq_name = T) + geom_seqlogo() + geom_msaBar()
pdf(file.path(base.dir, "msa", paste0("msa_","clustalW.pdf")))
temp
dev.off()
#subsetted regions 14150-16090
msa_default_align <- msaConvert(msa_default, type="seqinr::alignment")
temp <-ggmsa(msa_default_align,start=14150, end=16090, seq_name = T) #+ geom_seqlogo() + geom_msaBar()
pdf(file.path(base.dir, "msa", paste0("msa_","clustalW_14150-16090.pdf")))
temp
dev.off()
temp <-ggmsa(msa_default@unmasked,start=12920, end=18260, seq_name = T) + geom_seqlogo() + geom_msaBar()
pdf(file.path(base.dir, "msa", paste0("msa_","clustalW_12920-18260.pdf")))
temp
dev.off()

pdf(file.path(base.dir, "msa", "msa_ltr12_clustalW.pdf"))
msaPrettyPrint(msa_defaultmsa_default, output="asis", #y=c(164, 213),
    showNames="none", shadingMode="functional",
    shadingModeArg="structure",
    askForOverwrite=FALSE)
dev.off()

#same on ltr12 subsets
repeat_ltr12$repName <- droplevels(repeat_ltr12$repName )
repeat_ltr12_split<- split(repeat_ltr12, repeat_ltr12$repName)
ltr12_sequences_split <- lapply(repeat_ltr12_split, function(x)getSeq(genome,x))

#run msa with default sequences
dir.create(file.path(base.dir, "msa","msa_default_split"))
msa_default_split <- lapply(ltr12_sequences_split, function(x)msa(x, method="ClustalW", verbose=TRUE))
saveRDS(msa_default_split, file.path(base.dir, "msa","msa_default_split", "msa_default_split.rds"))
#msa_default_split <- readRDS(file.path(base.dir, "msa","msa_default_split", "msa_default_split.rds"))
#get consensus matrices
conMat_split <- lapply(msa_default_split, function(x)consensusMatrix(x))

#save alignments
for(i in names(msa_default_split)){
writeXStringSet(msa_default_split[[i]]@unmasked,
    file.path(base.dir, "msa","msa_default_split", paste0(i, "_msa_default_unmasked.fasta")))
}

#plot pretty msa
for(i in names(msa_default_split)){
#pdf(file.path(base.dir, "msa","msa_default_split", paste0("msa_",i,"_clustalW.pdf")))
    #msaPrettyPrint(msa_default_split[[i]],
    # showNames="none", shadingMode="functional",
    # shadingModeArg="structure",file=file.path(base.dir, "msa","msa_default_split", paste0("msa_",i,"_clustalW.asis")), output="asis", 
    # askForOverwrite=FALSE)
    #tools::texi2pdf(file.path(base.dir, "msa","msa_default_split", paste0("msa_",i,"_clustalW.tex")), clean=TRUE)
    #temp   <- msaConvert(msa_default_split[[i]], "bios2mds::align")
    #bios2mds::export.fasta(temp, outfile =file.path(base.dir, "msa","msa_default_split", paste0("msa_",i,"_clustalW.fa")), ncol = 60, open = "w")
    pdf(file.path(base.dir, "msa","msa_default_split", paste0("msa_",i,"_clustalW.pdf")))
        print(ggmsa(msa_default_split[[i]]@unmasked, seq_name = T) + geom_seqlogo() + geom_msaBar())
    dev.off()
    print(i)
    #dev.off()
}
    pdf(file.path(base.dir, "msa","msa_default_split", paste0("msa_",i,"_clustalW.pdf")))
ggmsa( msa_default_split[[i]]@unmasked, seq_name = T0) + geom_seqlogo() + geom_msaBar()
    dev.off()
