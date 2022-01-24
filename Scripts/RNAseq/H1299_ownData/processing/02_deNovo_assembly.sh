
#instal matching stringtie environment (same as long read assembly)
#conda create --name stringtie.2.1.1
conda activate stringtie.2.1.1
#conda install stringtie==2.1.1

# De novo assembly with stringtie
cd /omics/groups/OE0219/internal/tinat/210305_RNAseq_processing_LTR12Ctinat
mkdir /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4
mkdir /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/assembly/

#run de novo assembly for each bam file
conda activate stringtie==2.1.1
cd HISAT2/aligned_sorted/
for file in `ls *.bam`
do stringtie  ${file} --rf -o /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/assembly/${file}.gtf -p 5 \
-m 200 -f 0.05 -c 1.5 -p 10
echo $file
done

#merge de novo assembled trancripts
cd /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4
stringtie ./assembly/*.gtf --merge -G /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-p 10 -o ./mergedTranscripts.gtf 

#gffcompare 
conda activate gffcompare
mkdir /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffcompare
gffcompare -R -r /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-o /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/mergedTranscripts.gtf   
#sort gff compare
/home/heyj/tools/gff3sort/gff3sort.pl /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare/gffCompare.annotated.gtf  >/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare/gffCompare.annotated.sorted.gtf


#run ORF prediction
base.dir<- "/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef_analysis/"
#extract sequences of transcripts
gffread -w /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.fasta \
> -g /omics/groups/OE0219/internal/tinat/raw_data_repo/references/hg19.fa /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf
#get orfs
#final: top strand and 8 aa
cd tools/TransDecoder/
./TransDecoder.LongOrfs -t  /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.fasta \
> -m 8 -S -O /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa