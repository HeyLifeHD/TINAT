
#instal matching stringtie environment (same as long read assembly)
#conda create --name stringtie.2.1.1
conda activate stringtie.2.1.1
#conda install stringtie==2.1.1

# De novo assembly with stringtie
cd /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_knownRef/ 
mkdir -p /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/assembly/

#run de novo assembly for each bam file
cd HISAT2/aligned_sorted/
for file in `ls *.bam`
do stringtie  ${file} --rf -o /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/assembly/${file}.gtf \
-m 200 -f 0.01 -c 1 -p 15
echo $file
done

#merge de novo assembled trancripts
cd /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo
#for cellline in `ls ./assembly | tr '_' '\t'  | cut -f1 | sort | uniq`;do 
#echo ./assembly/${cellline}*.gtf
stringtie ./assembly/*.gtf --merge -G /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-p 15 -o ./SRP276887_mergedTranscripts.gtf 
#done

#gffcompare 
conda activate gffcompare
#for cellline in `ls ./assembly | tr '_' '\t'  | cut -f1 | sort | uniq`;do 
gffcompare -R -r /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-o /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare_SRP276887_mergedTranscripts /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/SRP276887_mergedTranscripts.gtf   
#echo $cellline 
#done

#sort gff compare
#for cellline in `ls ./assembly | tr '_' '\t'  | cut -f1 | sort | uniq`;do 
/home/heyj/tools/gff3sort/gff3sort.pl /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare_SRP276887_mergedTranscripts.annotated.gtf  > /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare_SRP276887_mergedTranscripts.annotated.sorted.gtf
#echo $cellline 
#done











#run ORF prediction
base.dir<- "/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef_analysis/"


#extract sequences of transcripts
gffread -w /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare.annotated.sorted.fasta \
> -g /omics/groups/OE0219/internal/tinat/raw_data_repo/references/hg19.fa /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare.annotated.sorted.gtf
#get orfs
#final: top strand and 8 aa
cd tools/TransDecoder/
./TransDecoder.LongOrfs -t  /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/gffCompare.annotated.sorted.fasta \
> -m 8 -S -O /omics/groups/OE0219/internal/tinat/220111_SRP276887_processing_deNovo/transdecoder_default_topStrand_8aa