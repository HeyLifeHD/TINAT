#Download  GSE140610 data (LPS stimulation of macrophages)
conda activate pysradb
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP276887
pysradb download -p SRP276887 --out-dir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP276887/
pysradb metadata SRP276887 > /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP276887/sample_sheet.txt

for file in `ls /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP276887/SRP276887/SRX109180{86..93}/*.sra`; do
numLines=$(fastq-dump -X 1 -Z --split-spot $file | wc -l)
#echo $file
if [ $numLines -eq 4 ]
then
    echo "$file is single-end"/
    fastq-dump ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP276887/single_end/
else
    echo "$file is paired-end"
    fastq-dump --split-files ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP276887/paired_end/
fi
done

cd /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP276887/paired_end/
parallel --gnu -j-2 --eta gzip ::: *.fastq
