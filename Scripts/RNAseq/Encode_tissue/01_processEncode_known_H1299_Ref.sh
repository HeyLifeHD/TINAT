module load anaconda3/2019.07
#create new nextflow environment
source activate nextflow_v21
#run this command a single time nextflow self-update
mkdir  /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing
cd  /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing
#copy or create rnaseq.config
cp  /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing_test/rnaseq.config /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing/
#create  design file
# export NXF_HOME=/omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing
# nextflow run nf-core/rnaseq -profile singularity -c rnaseq.config \
#     --input /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing/design.csv \
#     --outdir /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing/ \
#     -resume --skip_genebody_coverage


#same for H1299 reference genome
#run this command a single time nextflow self-update
mkdir  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref
cd  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref
#copy or create rnaseq.config
cp  /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing_test/rnaseq.config /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref/
cp  /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing/design.csv /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref/

export NXF_HOME=/omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref
nextflow run nf-core/rnaseq -profile singularity -c rnaseq.config \
    --input /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref/design.csv \
    --outdir /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref/ \
    --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
    --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
    -resume --skip_genebody_coverage --aligner hisat2 
    

#try with older version
mkdir  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon
cd  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon
#copy or create rnaseq.config
cp  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref/rnaseq.config /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/
nano /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/rnaseq.config
cp  /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing/design.csv /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/

export NXF_HOME=/omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon
nextflow run nf-core/rnaseq -profile singularity -c rnaseq.config \
    --input /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/design.csv \
    --outdir /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/ \
    --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
    --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
    -resume --skip_genebody_coverage --aligner hisat2 


#try with star salmon
#try with older version
mkdir  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon
cd  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon
#copy or create rnaseq.config
cp  /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref/rnaseq.config /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/
nano /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/rnaseq.config
cp  /omics/groups/OE0219/internal/tinat/Encode/220114_Encode_processing/design.csv /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/

export NXF_HOME=/omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon
nextflow run nf-core/rnaseq -profile singularity -c rnaseq.config \
    --input /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/design.csv \
    --outdir /omics/groups/OE0219/internal/tinat/Encode/220117_Encode_processing_H1299Ref_starSalmon/ \
    --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
    --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
    -resume --skip_genebody_coverage --aligner star_salmon 