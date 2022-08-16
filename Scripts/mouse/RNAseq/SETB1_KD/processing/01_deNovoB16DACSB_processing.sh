
#run processing
mkdir -p /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing
cd /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing

#copy data
mkdir /omics/groups/OE0219/internal/tinat/mouse_project/data/RNAseqSETB1KD_B16/
cp c010-datasets/Internal/Joschka/SetDB1/*/Raw*/*gz /omics/groups/OE0219/internal/tinat/mouse_project/data/RNAseqSETB1KD_B16/
#load anaconda module
module load anaconda3/2019.07
#load pipelin environment
source activate /omics/groups/OE0219/internal/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2

#copy and edit config file
cp /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_proecessing/project_specific.config \
/omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/project_specific.config
nano /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/project_specific.config

#run pipeline
/omics/groups/OE0219/internal/nextflow/nextflow \
-c /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/project_specific.config \
run /omics/groups/OE0219/internal/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing \
--reverseStranded  --pairedEnds --fasta /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/hisat2/genome.fa \
--gtf /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.gtf \
--hisat2_index /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/hisat2/ \
-profile cluster_odcf --skip_rseqc -resume --skip_genebody_coverage 


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
mkdir  /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/bigwig/
for file in `ls  /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/HISAT2/aligned_sorted/*.bam`;do
    pathname="${file%.*}"
    name=`basename $pathname`
    bamCoverage -b $file \
        -o  /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/bigwig/${name}.normalized.bw \
        -p 15 --normalizeUsing RPKM
    echo $name
done

#merge bigwigs
conda deactivate
conda activate kentUtils
cd  /omics/groups/OE0219/internal/tinat/mouse_project/220810_RNAseqSETB1KD_deNovoB16DACSB_proecessing/bigwig/

#control
bigWigMerge B16_control_sg1_rep1_1.sorted.normalized.bw B16_control_sg1_rep2_1.sorted.normalized.bw B16_control_sg2_rep1_1.sorted.normalized.bw B16_control_sg2_rep2_1.sorted.normalized.bw control.bedGraph
bedSort control.bedGraph control.bedGraph
bedGraphToBigWig control.bedGraph \
/omics/groups/OE0219/internal/genomes/Mmusculus/mm10/seq/mm10.chrom.sizes \
     control.bigWig 

#SETB1KD
bigWigMerge B16_SETDB1_KO_sg1_rep1_1.sorted.normalized.bw B16_SETDB1_KO_sg1_rep2_1.sorted.normalized.bw B16_SETDB1_KO_sg2_rep1_1.sorted.normalized.bw B16_SETDB1_KO_sg2_rep2_1.sorted.normalized.bw SETDB1_KO.bedGraph
bedSort SETDB1_KO.bedGraph SETDB1_KO.bedGraph
bedGraphToBigWig SETDB1_KO.bedGraph \
/omics/groups/OE0219/internal/genomes/Mmusculus/mm10/seq/mm10.chrom.sizes \
     SETDB1_KO.bigWig 