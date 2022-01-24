#run canonical processing of short read data with nf-core pipeline
module load anaconda3/2019.07
#load pipelin environment
source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
mkdir /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe
cd /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe


cp /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/project_specific.config \
 ./project_specific.config
nano ./project_specific.config

#start pipeline
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c  /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/ \
--fasta /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.fa --pairedEnd \
--gtf /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare_SRP276887_mergedTranscripts.annotated.sorted.gtf \
-profile cluster_odcf  --skip_genebody_coverage # -resume


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


#create bigwigs
conda activate deeptools
mkdir /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/bigwig/
for file in `ls /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/HISAT2/aligned_sorted/*.bam`;do
    pathname="${file%.*}"
    name=`basename $pathname`
    bamCoverage -b $file \
        -o /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/bigwig/${name}.bw \
        -p 10
    echo $name
done

#merge bigwigs
conda deactivate
conda activate kentUtils
cd /omics/groups/OE0219/internal/tinat/220111_SRP276887_deNovo_quantification_pe/bigwig/

#SETB1KO
bigWigMerge SRR1457491{6..9}_1.sorted.bw SETB1KO.bedGraph
bedSort SETB1KO.bedGraph SETB1KO.bedGraph
bedGraphToBigWig SETB1KO.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     SETB1KO.bigWig

#WT
bigWigMerge SRR1457492{0..3}_1.sorted.bw WT.bedGraph
bedSort WT.bedGraph WT.bedGraph
bedGraphToBigWig WT.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     WT.bigWig

# #split in forwad and reverse
# conda activate deeptools
# mkdir /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/
# for file in `ls /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/HISAT2/aligned_sorted/*.bam`;do
#     pathname="${file%.*}"
#     name=`basename $pathname`
#     bamCoverage -b $file \
#         -o /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/${name}_fw.bw \
#         -p 15 --filterRNAstrand forward --smoothLength 30 --binSize 10
#     bamCoverage -b $file \
#         -o /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/${name}_rv.bw \
#         -p 15 --filterRNAstrand reverse --smoothLength 30 --binSize 10
#     echo $name
# done

# #merge for forward
# #merge bigwigs
# conda deactivate
# conda activate kentUtils
# cd /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/

# #DMSO
# bigWigMerge AS-641773-LR-57123.sorted_fw.bw AS-641781-LR-57123.sorted_fw.bw AS-641789-LR-57123.sorted_fw.bw DMSO_fw.bedGraph
# bedSort DMSO_fw.bedGraph DMSO_fw.bedGraph
# bedGraphToBigWig DMSO_fw.bedGraph \
#  /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
#      DMSO_fw.bigWig