
#run processing
mkdir /omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/ 
cd /omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/ 
#load anaconda module
module load anaconda3/2019.07
#load pipelin environment
source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
cp  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/project_specific.config ./project_specific.config
nano  ./project_specific.config

#run pipeline
cd  /omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/ 
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c  /omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/  \
--reverseStranded  --pairedEnds --genome Hg19 \
-profile cluster_odcf  --skip_genebody_coverage #--pseudo aligner salmon

#copy transcript files into their own folder named by sample
cd ./stringtieFPKM/transcripts
mkdir ./../../transcripts
for file in *.gtf; do
    filename="${file%.*}"
    echo $filename
    mkdir  ../../transcripts/$filename
    cp $file ../../transcripts/$filename/
done
rmdir ./../../transcripts/ls

#run python script to create counts
cd ./../../transcripts
/omics/groups/OE0219/internal//nextflow/prepDE.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 
