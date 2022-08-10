#run canonical processing of short read data with nf-core pipeline
module load anaconda3/2019.07
#load pipelin environment
source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
mkdir /omics/groups/OE0219/internal/tinat/220118_SRP276887_deNovo_quantification_H12Ref
cd /omics/groups/OE0219/internal/tinat/220118_SRP276887_deNovo_quantification_H12Ref


cp /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/project_specific.config \
 ./project_specific.config
nano ./project_specific.config

#start pipeline
/omics/groups/OE0219/internal/nextflow/nextflow \
    -c  /omics/groups/OE0219/internal/tinat/220118_SRP276887_deNovo_quantification_H12Ref/project_specific.config \
    run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
    --outdir /omics/groups/OE0219/internal/tinat/220118_SRP276887_deNovo_quantification_H12Ref/ \
    --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa --pairedEnd \
    --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
    -profile cluster_odcf  --skip_genebody_coverage # -resume


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
/omics/groups/OE0219/internal/nextflow/prepDE.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 

