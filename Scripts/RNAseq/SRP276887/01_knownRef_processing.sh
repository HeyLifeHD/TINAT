
#run processing
mkdir /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/ 
cd /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/ 
#load anaconda module
module load anaconda3/2019.07
#load pipelin environment
source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
cp /omics/groups/OE0219/internal/repguide/201207_RNAseq_processing_hg19_incl_gRNAseq/project_specific.config \
/omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/project_specific.config
nano /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/project_specific.config

#run pipeline
cd /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/ \
--pairedEnd --genome Hg19 \
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

#run python script to create counts
cd ./../../transcripts+
/omics/groups/OE0219/internal/nextflow/prepDE.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 120 
