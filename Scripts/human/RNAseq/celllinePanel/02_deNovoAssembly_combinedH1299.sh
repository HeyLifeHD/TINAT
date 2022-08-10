
#copy h1299 and novel cellline data into common folder
mkdir -p /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/HISAT2/aligned_sorted/
cp /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification/HISAT2/aligned_sorted/A* \
/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/HISAT2/aligned_sorted/
##copy missing
cp /omics/groups/OE0219/internal/tinat/Cellline_panel/220309_deNovo_quantification_H1299ref_reverse/HISAT2/aligned_sorted/A* \
/omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/HISAT2/aligned_sorted/


# De novo assembly with stringtie
cd /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/
mkdir ./assembly/

#run de novo assembly for each bam file
conda activate stringtie.2.1.1
cd HISAT2/aligned_sorted/
for file in `ls *.bam`
do stringtie  ${file} --rf -o /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/assembly/${file}.gtf \
-m 200 -f 0.05 -c 1.5 -p 10
echo $file
done

#merge de novo assembled trancripts
cd /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/
stringtie ./assembly/*.gtf --merge -G /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-p 10 -o ./mergedTranscripts.gtf 

#gffcompare 
conda activate gffcompare
gffcompare -R -r /omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
-o /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/mergedTranscripts.gtf   
#sort gff compare
/home/heyj/tools/gff3sort/gff3sort.pl /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.gtf  > /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.gtf

#subset based on strnads
awk '($7 == "+" || $7 == "-")' /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.gtf > /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted_sub.gtf

#run ORF prediction
#extract sequences of transcripts
gffread -w /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.fasta \
-g /omics/groups/OE0219/internal/tinat/raw_data_repo/references/hg19.fa /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.gtf
#get orfs
#final: top strand and 8 aa
cd /home/heyj/tools/TransDecoder/
./TransDecoder.LongOrfs -t  /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.fasta \
-m 8 -S -O /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/transdecoder_default_topStrand_8aa
