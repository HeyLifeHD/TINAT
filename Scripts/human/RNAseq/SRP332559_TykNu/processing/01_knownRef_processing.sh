
#run processing
mkdir /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ 
cd /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ 
#load anaconda module
module load anaconda3/2019.07
#load pipelin environment
source activate /icgc/dkfzlsdf/analysis/C010/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
cp /omics/groups/OE0219/internal/tinat/201207_RNAseq_processing_hg19_incl_gRNAseq/project_specific.config \
 /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/project_specific.config
nano  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /project_specific.config
#run pipeline
cd  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ 
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow \
-c  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /project_specific.config \
run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/  \
--reverseStranded  --pairedEnds --genome Hg19 \
-profile cluster_odcf  --skip_genebody_coverage  -resume #- #--pseudo aligner salmon
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


#create bigwigs
conda activate deeptools
mkdir  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/
for file in `ls  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /HISAT2/aligned_sorted/*.bam`;do
    pathname="${file%.*}"
    name=`basename $pathname`
    bamCoverage -b $file \
        -o  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/${name}.bw \
        -p 15
    echo $name
done

#merge bigwigs
conda deactivate
conda activate kentUtils
cd  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/

#DMSO
bigWigMerge AS-641773-LR-57123.sorted.bw AS-641781-LR-57123.sorted.bw AS-641789-LR-57123.sorted.bw DMSO.bedGraph
bedSort DMSO.bedGraph DMSO.bedGraph
bedGraphToBigWig DMSO.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DMSO.bigWig

#DAC
bigWigMerge AS-641775-LR-57123.sorted.bw AS-641783-LR-57123.sorted.bw AS-641791-LR-57123.bw DAC.bedGraph
bedSort DAC.bedGraph DAC.bedGraph
bedGraphToBigWig DAC.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DAC.bigWig

#SB939
bigWigMerge AS-641777-LR-57123.sorted.bw AS-641785-LR-57123.bw AS-641793-LR-57123.sorted.bw SB939.bedGraph
bedSort SB939.bedGraph SB939.bedGraph
bedGraphToBigWig LTR12C_wG.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     SB939.bigWig


#DACSB
bigWigMerge AS-641779-LR-57123.sorted.bw AS-641787-LR-57123.bw AS-641795-LR-57123.sorted.bw DACSB939.bedGraph
bedSort SB939.bedGraph DACSB939.bedGraph
bedGraphToBigWig DACSB939.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DACSB939.bigWig




#split in forwad and reverse
conda activate deeptools
mkdir  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/
for file in `ls  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /HISAT2/aligned_sorted/*.bam`;do
    pathname="${file%.*}"
    name=`basename $pathname`
    bamCoverage -b $file \
        -o  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/${name}_fw.bw \
        -p 15 --filterRNAstrand forward --smoothLength 30 --binSize 10
    bamCoverage -b $file \
        -o  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/${name}_rv.bw \
        -p 15 --filterRNAstrand reverse --smoothLength 30 --binSize 10
    echo $name
done

#merge for forward
#merge bigwigs
conda deactivate
conda activate kentUtils
cd  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/

#DMSO
bigWigMerge AS-641773-LR-57123.sorted_fw.bw AS-641781-LR-57123.sorted_fw.bw AS-641789-LR-57123.sorted_fw.bw DMSO_fw.bedGraph
bedSort DMSO_fw.bedGraph DMSO_fw.bedGraph
bedGraphToBigWig DMSO_fw.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DMSO_fw.bigWig

#DAC
bigWigMerge AS-641775-LR-57123.sorted_fw.bw AS-641783-LR-57123.sorted_fw.bw AS-641791-LR-57123_fw.bw DAC_fw.bedGraph
bedSort DAC_fw.bedGraph DAC_fw.bedGraph
bedGraphToBigWig DAC_fw.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DAC_fw.bigWig

#SB939
bigWigMerge AS-641777-LR-57123_fw.sorted.bw AS-641785-LR-57123_fw.bw AS-641793-LR-57123.sorted_fw.bw SB939_fw.bedGraph
bedSort SB939_fw.bedGraph SB939_fw.bedGraph
bedGraphToBigWig LTR12C_wG_fw.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     SB939_fw.bigWig


#DACSB
bigWigMerge AS-641779-LR-57123_fw.sorted.bw AS-641787-LR-57123_fw.bw AS-641795-LR-57123_fw.sorted.bw DACSB939_fw.bedGraph
bedSort SB939_fw.bedGraph DACSB939_fw.bedGraph
bedGraphToBigWig DACSB939.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DACSB939_fw.bigWig



#merge for reverse
#merge bigwigs
conda deactivate
conda activate kentUtils
cd  /omics/groups/OE0219/internal/tinat/210920_SRP332559_processing_knownRef/ /bigwig/

#DMSO
bigWigMerge AS-641773-LR-57123.sorted_rv.bw AS-641781-LR-57123.sorted_rv.bw AS-641789-LR-57123.sorted_rv.bw DMSO_rv.bedGraph
bedSort DMSO_rv.bedGraph DMSO_rv.bedGraph
bedGraphToBigWig DMSO_rv.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DMSO_rv.bigWig

#DAC
bigWigMerge AS-641775-LR-57123.sorted_rv.bw AS-641783-LR-57123.sorted_rv.bw AS-641791-LR-57123_rv.bw DAC_rv.bedGraph
bedSort DAC_rv.bedGraph DAC_rv.bedGraph
bedGraphToBigWig DAC_rv.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DAC_rv.bigWig

#SB939
bigWigMerge AS-641777-LR-57123_rv.sorted.bw AS-641785-LR-57123_rv.bw AS-641793-LR-57123.sorted_rv.bw SB939_rv.bedGraph
bedSort SB939_rv.bedGraph SB939_rv.bedGraph
bedGraphToBigWig LTR12C_wG_rv.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     SB939_rv.bigWig


#DACSB
bigWigMerge AS-641779-LR-57123_rv.sorted.bw AS-641787-LR-57123_rv.bw AS-641795-LR-57123_rv.sorted.bw DACSB939_rv.bedGraph
bedSort SB939_rv.bedGraph DACSB939_rv.bedGraph
bedGraphToBigWig DACSB939.bedGraph \
 /icgc/dkfzlsdf/analysis/C010/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DACSB939_rv.bigWig