
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




#get normalized bigwigs
#create bigwigs
conda activate deeptools
mkdir  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/
for file in `ls /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/HISAT2/aligned_sorted/*.bam`;do
    pathname="${file%.*}"
    name=`basename $pathname`
    bamCoverage -b $file \
        -o  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/${name}.normalized.bw \
        -p 15 --normalizeUsing RPKM
    echo $name
done


#merge bigwigs
conda deactivate
conda activate kentUtils
cd  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/bigwig/

#DMSO
bigWigMerge -adjust=0.3333 AS-811002-LR-62189.sorted.normalized.bw AS-811010-LR-62189.sorted.normalized.bw AS-811018-LR-62189.sorted.normalized.bw DMSO.normalized.bedGraph
bedSort DMSO.normalized.bedGraph DMSO.normalized.bedGraph
bedGraphToBigWig DMSO.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/seq/mm10.chrom.sizes \
     DMSO.normalized.bigWig

#DAC
bigWigMerge -adjust=0.5 AS-811012-LR-62189.sorted.normalized.bw AS-811020-LR-62189.sorted.normalized.bw DAC.normalized.bedGraph
bedSort DAC.normalized.bedGraph DAC.normalized.bedGraph
bedGraphToBigWig DAC.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/seq/mm10.chrom.sizes \
     DAC.normalized.bigWig

#SB939
bigWigMerge -adjust=0.3333 AS-811006-LR-62189.sorted.normalized.bw AS-811014-LR-62189.sorted.normalized.bw AS-811022-LR-62189.sorted.normalized.bw SB939.normalized.bedGraph
bedSort SB939.normalized.bedGraph SB939.normalized.bedGraph
bedGraphToBigWig SB939.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/seq/mm10.chrom.sizes \
     SB939.normalized.bigWig


#DACSB
bigWigMerge  -adjust=0.5 AS-811008-LR-62189.sorted.normalized.bw AS-811016-LR-62189.sorted.normalized.bw  DACSB939.normalized.bedGraph
bedSort DACSB939.normalized.bedGraph DACSB939.normalized.bedGraph
bedGraphToBigWig DACSB939.normalized.bedGraph \
 /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/seq/mm10.chrom.sizes \
     DACSB939.normalized.bigWig

