
#instal matching stringtie environment (same as long read assembly)
#conda create --name stringtie.2.1.1
conda activate stringtie.2.1.1
#conda install stringtie==2.1.1

# De novo assembly with stringtie
cd /omics/groups/OE0219/internal/tinat/mouse_project/220808_RNAseq_knownGene_proecessing
mkdir -p /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/assembly/

#run de novo assembly for each bam file
conda activate stringtie==2.1.1
cd HISAT2/aligned_sorted/

#remove AS-811024-LR-62189 from analysis
mkdir failed
mv AS-811024-LR-62189* failed/
for file in `ls *.bam`
do stringtie  ${file} --rf -o /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/assembly/${file}.gtf \
-m 200 -f 0.05 -c 1.5 -p 10
echo $file
done

#merge de novo assembled trancripts
cd /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo
stringtie ./assembly/*.gtf --merge -G /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/rnaseq/gencode.vM19.annotation.gtf \
-p 10 -o ./mergedTranscripts.gtf 

#gffcompare 
conda activate gffcompare
gffcompare -R -r /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/rnaseq/gencode.vM19.annotation.gtf \
-o /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/mergedTranscripts.gtf   
#sort gff compare
/home/heyj/tools/gff3sort/gff3sort.pl /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.gtf  > /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.gtf

#subset based on strnads
#/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.gtf
awk '($7 == "+" || $7 == "-")' /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.gtf > /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted_sub.gtf










#run ORF prediction
#extract sequences of transcripts
gffread -w /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.fasta \
> -g /omics/groups/OE0219/internal/genomes/Mmusculus/mm10/hisat2/genome.fa /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.gtf
#get orfs
#final: top strand and 8 aa
cd tools/TransDecoder/
./TransDecoder.LongOrfs -t  /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.fasta \
> -m 8 -S -O /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/transdecoder_default_topStrand_8aa_test


mkdir /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/transdecoder_default_topStrand_8aa_test
cd /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/transdecoder_default_topStrand_8aa_test
/home/heyj/tools/TransDecoder/TransDecoder.LongOrfs -t /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted.fasta -m 8 -S -O /omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/transdecoder_default_topStrand_8aa_test/