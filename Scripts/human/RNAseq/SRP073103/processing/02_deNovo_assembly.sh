
#instal matching stringtie environment (same as long read assembly)
#conda create --name stringtie.2.1.1
conda activate stringtie.2.1.1
#conda install stringtie==2.1.1

# De novo assembly with stringtie
cd /omics/groups/OE0219/internal/tinat/210913_SRP073103_processing_knownRef/
mkdir -p /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines
mkdir /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/assembly/

#run de novo assembly for each bam file
conda activate stringtie==2.1.1
cd HISAT2/aligned_sorted/
for file in `ls *.bam`
do stringtie  ${file} --rf -o /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/assembly/${file}.gtf -p 5 \
-m 200 -f 0.01 -c 1 -p 10
echo $file
done

#merge de novo assembled trancripts
cd /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines
for cellline in `ls ./assembly | tr '_' '\t'  | cut -f1 | sort | uniq`;do 
ls  ./assembly/${cellline}*.gtf
stringtie ./assembly/${cellline}*.gtf --merge -G /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-p 10 -o ./${cellline}_mergedTranscripts.gtf 
done

#gffcompare 
conda activate gffcompare
for cellline in `ls ./assembly | tr '_' '\t'  | cut -f1 | sort | uniq`;do 
gffcompare -R -r /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-o /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/gffCompare_${cellline} /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/${cellline}_mergedTranscripts.gtf   
echo $cellline 
done

#sort gff compare
for cellline in `ls ./assembly | tr '_' '\t'  | cut -f1 | sort | uniq`;do 
/home/heyj/tools/gff3sort/gff3sort.pl /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/gffCompare_${cellline}.annotated.gtf  > /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/gffCompare_${cellline}.annotated.sorted.gtf
echo $cellline 
done


#run ORF prediction
base.dir<- "/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef_analysis/"


#extract sequences of transcripts
gffread -w /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/gffCompare.annotated.sorted.fasta \
> -g /omics/groups/OE0219/internal/tinat/raw_data_repo/references/hg19.fa /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/gffCompare.annotated.sorted.gtf
#get orfs
#final: top strand and 8 aa
cd tools/TransDecoder/
./TransDecoder.LongOrfs -t  /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/gffCompare.annotated.sorted.fasta \
> -m 8 -S -O /omics/groups/OE0219/internal/tinat//210913_SRP073103_processing_deNovo/all_celllines/transdecoder_default_topStrand_8aa