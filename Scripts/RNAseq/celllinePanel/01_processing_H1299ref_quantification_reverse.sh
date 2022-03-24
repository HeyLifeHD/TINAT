#run canonical processing of short read data with nf-core pipeline
module load anaconda3/2019.07
#load pipelin environment
source activate /omics/groups/OE0219/internal/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
mkdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse
cd /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification
cp /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification/project_specific.config  \
 ./project_specific.config
nano ./project_specific.config

#start pipeline
/omics/groups/OE0219/internal/nextflow/nextflow -c /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse/project_specific.config \
run /omics/groups/OE0219/internal/nextflow/rnaseq-master/main.nf \
--outdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse/ \
--reverseStranded --pairedEnds  \
--fasta                                 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
--gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
--hisat2_index /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification/hisat2_index/ \
-resume -profile cluster_odcf  #--skip_genebody_coverage 

cd /omics/groups/OE0219/internal/tinat/220309_deNovo_quantification_H1299ref_reverse_noDup
/omics/groups/OE0219/internal/nextflow/nextflow -c /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse/project_specific.config \
run /omics/groups/OE0219/internal/nextflow/rnaseq-master/main.nf \
--outdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse_noDup/ \
--reverseStranded --pairedEnds  \
--fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
--gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
--hisat2_index /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification/hisat2_index/ \
 --skip_markduplicates -profile cluster_odcf  #--skip_genebody_coverage 




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




#same for newer pipeline version
mkdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse_v14
cd /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse_v14
cp  /omics/groups/OE0219/internal/tinat/Encode/220127_deNovo_quantification_H1299ref_reverse/rnaseq.config ./
#nano rnaseq.config 

export NXF_HOME=/omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse_v14
nextflow run nf-core/rnaseq \
    -resume -profile singularity \
    -c rnaseq.config \
    -r 1.4 \
    --aligner hisat2 \
    --reads '/omics/odcf/project/OE0219/tinat/sequencing/rna_sequencing/core/run22*/*{1,2}.fastq.gz' \
    --outdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse_v14/ \
    --pairedEnd \
    --reverse_stranded \
    --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
    --gtf /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf \
    --email j.hey@dkfz.de 

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



