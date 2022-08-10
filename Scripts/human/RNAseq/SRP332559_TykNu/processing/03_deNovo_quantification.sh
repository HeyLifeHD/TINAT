#run canonical processing of short read data with nf-core pipeline
module load anaconda3/2019.07
#load pipelin environment
source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
mkdir /omics/groups/OE0219/internal/tinat/210921_SRP332559_deNovo_quantification
cd /omics/groups/OE0219/internal/tinat/210921_SRP332559_deNovo_quantification

cp /icgc/dkfzlsdf/analysis/C010/repguide/201207_RNAseq_processing_hg19_incl_gRNAseq/project_specific.config \
 ./project_specific.config
 nano ./project_specific.config

#start pipeline
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c  /omics/groups/OE0219/internal/tinat/210921_SRP332559_deNovo_quantification/project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/210921_SRP332559_deNovo_quantification/ \
--fasta /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.fa \
--gtf /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_deNovo/gffCompare_SRP332559_mergedTranscripts.annotated.sorted.gtf \
--reverseStranded  --pairedEnds \
-profile cluster_odcf  --skip_genebody_coverage -resume


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
/omics/groups/OE0219/internal/nextflow/prepDE.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 

