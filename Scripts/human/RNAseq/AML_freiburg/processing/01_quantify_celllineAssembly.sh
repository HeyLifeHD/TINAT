#process rnaseq of aml patients from freibur
#copy data
cp /mnt/external_freiburg/RNAseq_patient_blasts/* /omics/groups/OE0219/internal/tinat/external_data/RNAseq_Freiburg/data/   
cp /mnt/external_freiburg/RNAseq_patient_blasts/'Galaxy10-[GG_A01032_c1d1_R1.fastq].fastqsanger' /omics/groups/OE0219/internal/tinat/external_data/RNAseq_Freiburg/data/   
#rename and pack data
cd /omics/groups/OE0219/internal/tinat/external_data/RNAseq_Freiburg/data/
for file in *; do
NAME=`echo "$file" | cut -d'[' -f2`
echo $NAME
NAMENEW=`echo "$NAME" | cut -d']' -f1`
echo $NAMENEW
mv $file $NAMENEW
done
gzip *




#run canonical processing of short read data with nf-core pipeline
module load anaconda3/2019.07
#load pipelin environment
source activate /omics/groups/OE0219/internal/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
mkdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220607_AMLPatient_deNovoCellline_quantification
cd /omics/groups/OE0219/internal/tinat/Cellline_panel/220607_AMLPatient_deNovoCellline_quantification/
cp /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification/project_specific.config  \
 ./project_specific.config
nano ./project_specific.config

#start pipeline
/omics/groups/OE0219/internal/nextflow/nextflow \
-c /omics/groups/OE0219/internal/tinat/Cellline_panel/220607_AMLPatient_deNovoCellline_quantification/project_specific.config \
run /omics/groups/OE0219/internal/nextflow/rnaseq-master/main.nf \
--outdir  /omics/groups/OE0219/internal/tinat/Cellline_panel/220607_AMLPatient_deNovoCellline_quantification/ \
--reverseStranded --pairedEnds  \
--fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
--gtf /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.gtf \
-profile cluster_odcf \
--hisat2_index /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification/hisat2_index/ \
-resume --skip_genebody_coverage --skip_markduplicates


# #singularity
# module load anaconda3/2019.07
# #create new nextflow environment
# source activate nextflow_v21
# #copy config
# mkdir  /omics/groups/OE0219/internal/tinat/Cellline_panel/220318_deNovo_quantification_H1299ref_reverse_quantification
# cd  /omics/groups/OE0219/internal/tinat/Cellline_panel/220318_deNovo_quantification_H1299ref_reverse_quantification
# nano ./rnaseq.config
# #nano rnaseq.config 

# #run piepeline
# export NXF_HOME=/omics/groups/OE0219/internal/tinat/Cellline_panel/220318_deNovo_quantification_H1299ref_reverse_quantification/
# nextflow run nf-core/rnaseq \
#     -profile singularity \
#     -c rnaseq.config \
#     -r 1.2 \
#     --aligner hisat2 \
#     --reads '/omics/odcf/project/OE0219/tinat/sequencing/rna_sequencing/core/run22*/*{1,2}.fastq.gz' \
#     --outdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220318_deNovo_quantification_H1299ref_reverse_quantification/ \
#     --pairedEnd \
#     --reverse_stranded \
#     --fasta /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.fa \
#     --gtf  /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.gtf \
#     #--skip_genebody_coverage -resume  \
#     --email j.hey@dkfz.de --max_memory 120.GB --max_time '48.h' --max_cpus 16




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



