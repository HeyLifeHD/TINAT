
#run processing
mkdir -p /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing
cd /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing

#load anaconda module
module load anaconda3/2019.07
#load pipelin environment
source activate /omics/groups/OE0219/internal/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2

#copy and edit config file
cp /omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing/project_specific.config \
/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/project_specific.config
nano /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/project_specific.config


#run pipeline
/omics/groups/OE0219/internal/nextflow/nextflow \
-c /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/project_specific.config \
run /omics/groups/OE0219/internal/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing \
--reverseStranded  --pairedEnds --fasta /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/hisat2/genome.fa \
--gtf /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.gtf \
--hisat2_index /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/hisat2/ \
-profile cluster_odcf  --skip_genebody_coverage --skip_rseqc -resume #--pseudo aligner salmon


#copy transcript files into their own folder named by sample
cd ./stringtieFPKM/transcripts
mkdir ./../../transcripts
for file in *.gtf; do
    filename="${file%.*}"
    echo $filename
    mkdir  ../../transcripts/$filename
    cp $file ../../transcripts/$filename/
done

#run python script to create counts
cd ./../../transcripts
/omics/groups/OE0219/internal/nextflow/prepDE_v2.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 



#create bigwigs
conda activate deeptools
mkdir  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/
for file in `ls  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/HISAT2/aligned_sorted/*.bam`;do
    pathname="${file%.*}"
    name=`basename $pathname`
    bamCoverage -b $file \
        -o  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/${name}.bw \
        -p 15
    echo $name
done

#merge bigwigs
conda deactivate
conda activate kentUtils
cd  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/

#DMSO
bigWigMerge AS-641773-LR-57123.sorted.bw AS-641781-LR-57123.sorted.bw AS-641789-LR-57123.sorted.bw DMSO.bedGraph
bedSort DMSO.bedGraph DMSO.bedGraph
bedGraphToBigWig DMSO.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DMSO.bigWig

#DAC
bigWigMerge AS-641775-LR-57123.sorted.bw AS-641783-LR-57123.sorted.bw AS-641791-LR-57123.bw DAC.bedGraph
bedSort DAC.bedGraph DAC.bedGraph
bedGraphToBigWig DAC.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DAC.bigWig

#SB939
bigWigMerge AS-641777-LR-57123.sorted.bw AS-641785-LR-57123.bw AS-641793-LR-57123.sorted.bw SB939.bedGraph
bedSort SB939.bedGraph SB939.bedGraph
bedGraphToBigWig LTR12C_wG.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     SB939.bigWig


#DACSB
bigWigMerge AS-641779-LR-57123.sorted.bw AS-641787-LR-57123.bw AS-641795-LR-57123.sorted.bw DACSB939.bedGraph
bedSort SB939.bedGraph DACSB939.bedGraph
bedGraphToBigWig DACSB939.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DACSB939.bigWig



#get normalized bigwigs
#create bigwigs
conda activate deeptools
mkdir  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/
for file in `ls  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/HISAT2/aligned_sorted/*.bam`;do
    pathname="${file%.*}"
    name=`basename $pathname`
    bamCoverage -b $file \
        -o  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/${name}.normalized.bw \
        -p 15 --normalizeUsing RPKM
    echo $name
done


-------
#merge bigwigs
conda deactivate
conda activate kentUtils
cd  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/

#DMSO
bigWigMerge AS-641773-LR-57123.sorted.normalized.bw AS-641781-LR-57123.sorted.normalized.bw AS-641789-LR-57123.sorted.normalized.bw DMSO.normalized.bedGraph
bedSort DMSO.normalized.bedGraph DMSO.normalized.bedGraph
bedGraphToBigWig DMSO.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DMSO.normalized.bigWig

#DAC
bigWigMerge AS-641775-LR-57123.sorted.normalized.bw AS-641783-LR-57123.sorted.normalized.bw AS-641791-LR-57123.sorted.normalized.bw DAC.normalized.bedGraph
bedSort DAC.normalized.bedGraph DAC.normalized.bedGraph
bedGraphToBigWig DAC.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DAC.normalized.bigWig

#SB939
bigWigMerge AS-641777-LR-57123.sorted.normalized.bw AS-641785-LR-57123.sorted.normalized.bw AS-641793-LR-57123.sorted.normalized.bw SB939.normalized.bedGraph
bedSort SB939.normalized.bedGraph SB939.normalized.bedGraph
bedGraphToBigWig SB939.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     SB939.normalized.bigWig


#DACSB
bigWigMerge AS-641779-LR-57123.sorted.normalized.bw AS-641787-LR-57123.sorted.normalized.bw AS-641795-LR-57123.sorted.normalized.bw DACSB939.normalized.bedGraph
bedSort DACSB939.normalized.bedGraph DACSB939.normalized.bedGraph
bedGraphToBigWig DACSB939.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DACSB939.normalized.bigWig





