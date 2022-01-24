#run canonical processing of short read data with nf-core pipeline
module load anaconda3/2019.07
#load pipelin environment
source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
mkdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification

#run quantification for each cellline seperately 
#LNT-229
mkdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/LNT-229
cd /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/LNT-229

cp /icgc/dkfzlsdf/analysis/C010/repguide/201207_RNAseq_processing_hg19_incl_gRNAseq/project_specific.config \
 ./project_specific.config
 nano ./project_specific.config

#start pipeline
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c  /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/LNT-229/project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/LNT-229/ \
--fasta /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.fa --singleEnd \
--gtf /omics/groups/OE0219/internal/tinat/210913_SRP073103_processing_deNovo/all_celllines/gffCompare_LNT-229.annotated.sorted.gtf \
-profile cluster_odcf  #--skip_genebody_coverage -resume


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


#######/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/merged_fastq/LNT*.fastq.gz
#T98G
mkdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/T98G
cd /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/T98G

cp /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/LNT-229/project_specific.config \
 ./project_specific.config
 nano ./project_specific.config

#start pipeline
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c  /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/T98G/project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/T98G/ \
 --singleEnd  \
--fasta /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.fa \
--gtf /omics/groups/OE0219/internal/tinat/210913_SRP073103_processing_deNovo/all_celllines/gffCompare_T98G.annotated.sorted.gtf \
-profile cluster_odcf  #--skip_genebody_coverage -resume


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



#######
#U-87
mkdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/U-87
cd /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/U-87

cp /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/LNT-229/project_specific.config \
 ./project_specific.config
 nano ./project_specific.config

#start pipeline
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c  /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/U-87/project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_deNovo_quantification/U-87/ \
 --singleEnd    \
--fasta /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.fa \
--gtf /omics/groups/OE0219/internal/tinat/210913_SRP073103_processing_deNovo/all_celllines/gffCompare_U-87.annotated.sorted.gtf \
-profile cluster_odcf # --skip_genebody_coverage -resume


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





