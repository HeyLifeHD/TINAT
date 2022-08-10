
#run processing
mkdir -p/omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing
cd /omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing

#copy sequencing data
mkdir -p /omics/groups/OE0219/internal/tinat/mouse_project/data/RNAseq
cp /omics/odcf/project/OE0219/tinat/sequencing/rna_sequencing/view-by-pid/OE0219_TINAT_B16-F10/*/paired/run220524_VH00693_45_AAACLTYHV/sequence/*.fastq.gz /omics/groups/OE0219/internal/tinat/mouse_project/data/RNAseq/ 
AS-811004-LR-62189_R1.fastq.gz
#load anaconda module
module load anaconda3/2019.07
#load pipelin environment
source activate /omics/groups/OE0219/internal/nextflow/environments/rnaseq-master.env
#source activate nf-core-rnaseq-1.2
#copy and edit config file
cp /omics/groups/OE0219/internal/repguide/210305_RNAseq_processing_LTR12CRepguide/project_specific.config \
/omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing/project_specific.config
nano /omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing/project_specific.config


#run pipeline
/omics/groups/OE0219/internal/nextflow/nextflow \
-c/omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing
/project_specific.config \
run /omics/groups/OE0219/internal/nextflow/rnaseq-master/main.nf  \
--outdir /omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing \
--reverseStranded  --pairedEnds --fasta /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/hisat2/genome.fa \
--gtf /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/rnaseq/gencode.vM19.annotation.gtf \
--hisat2_index /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/hisat2/ \
-profile cluster_odcf  --skip_genebody_coverage --skip_rseqc -resume #--pseudo aligner salmon


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









