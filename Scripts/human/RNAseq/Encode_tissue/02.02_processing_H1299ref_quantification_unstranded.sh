#try with singularity
module load anaconda3/2019.07
#create new nextflow environment
source activate nextflow_v21
#copy config

cp  /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_reverse/rnaseq.config ./
nano rnaseq.config 
#run piepeline
export NXF_HOME=/omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded
nextflow run nf-core/rnaseq \
    -profile singularity \
    -c rnaseq.config \
    -r 1.4 \
    --aligner hisat2 \
    --reads '/omics/groups/OE0219/internal/tinat/Encode/data_new/unstranded/*{1,2}.fastq.gz' \
    --outdir /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded/ \
    --pairedEnd \
    --unstranded \
    --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
    --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
    --skip_genebody_coverage \
    --email j.hey@dkfz.de -resume --skipPreseq --skipBiotypeQC --skipDupRadar --skipQualimap --skipRseQC 


#copy transcript files into their own folder named by sample
cd ./stringtieFPKM/transcripts
mkdir ./../../transcripts
for file in *.gtf; do
    filename="${file%.*}"
    echo $file
    echo $filename
    mkdir  ../../transcripts/$filename
    cp $file ../../transcripts/$filename/
done

#run python script to create counts
cd ./../../transcripts
/omics/groups/OE0219/internal/nextflow/prepDE_v2.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 

# #run canonical processing of short read data with nf-core pipeline
# module load anaconda3/2019.07
# #load pipelin environment
# source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
# #source activate nf-core-rnaseq-1.2
# #copy and edit config file
# mkdir /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded
# cd /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded


# cp /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/project_specific.config \
#  ./project_specific.config
# nano ./project_specific.config

# #start pipeline
# /omics/groups/OE0219/internal/nextflow/nextflow \
#     -resume -c  /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded/project_specific.config \
#     run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
#     --outdir /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded/ \
#     --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa --pairedEnd \
#     --unstranded \
#     --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
#     -profile cluster_odcf  --skip_genebody_coverage 

# #copy transcript files into their own folder named by sample
# cd ./stringtieFPKM/transcripts
# mkdir ./../../transcripts
# for file in *.gtf; do
#     filename="${file%.*}"
#     echo $file
#     echo $filename
#     mkdir  ../../transcripts/$filename
#     cp $file ../../transcripts/$filename/
# done

# #run python script to create counts
# cd ./../../transcripts
# /omics/groups/OE0219/internal/nextflow/prepDE.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 

# #try with singularity
# module load anaconda3/2019.07
# #create new nextflow environment
# source activate nextflow_v21
# #copy config
# mkdir /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded_v12
# cd /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded_v12
# cp  /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_reverse_v12/rnaseq.config ./
# #nano rnaseq.config 

# #run piepeline
# export NXF_HOME=/omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded_v12
# nextflow run nf-core/rnaseq \
#     -profile singularity \
#     -c rnaseq.config \
#     -r 1.2 \
#     --aligner hisat2 \
#     --reads '/omics/groups/OE0219/internal/tinat/Encode/data_new/unstranded/*{1,2}.fastq.gz' \
#     --outdir /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_unstranded_v12/ \
#     --pairedEnd \
#     --unstranded \
#     --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
#     --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
#     --skip_genebody_coverage \
#     --email j.hey@dkfz.de --skipPreseq --skipBiotypeQC --skipDupRadar --skipQualimap --skipRseQC # -resume 


# #copy transcript files into their own folder named by sample
# cd ./stringtieFPKM/transcripts
# mkdir ./../../transcripts
# for file in *.gtf; do
#     filename="${file%.*}"
#     echo $file
#     echo $filename
#     mkdir  ../../transcripts/$filename
#     cp $file ../../transcripts/$filename/
# done

# #run python script to create counts
# cd ./../../transcripts
# /omics/groups/OE0219/internal/nextflow/prepDE.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 
