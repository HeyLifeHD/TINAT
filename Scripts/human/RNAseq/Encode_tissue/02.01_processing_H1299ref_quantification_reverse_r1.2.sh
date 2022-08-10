#singularity
module load anaconda3/2019.07
#create new nextflow environment
source activate nextflow_v21
#copy config
cd /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_reverse
cp  /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded/rnaseq.config ./
#nano rnaseq.config 

#run piepeline
export NXF_HOME=/omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_reverse
nextflow run nf-core/rnaseq \
    -resume -profile singularity \
    -c rnaseq.config \
    -r 1.2 \
    --aligner hisat2 \
    --reads '/omics/groups/OE0219/internal/tinat/Encode/data_new/reverse/*{1,2}.fastq.gz' \
    --outdir /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_reverse/ \
    --pairedEnd \
    --reverse_stranded \
    --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
    --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
    --skip_genebody_coverage \
    --email j.hey@dkfz.de 
