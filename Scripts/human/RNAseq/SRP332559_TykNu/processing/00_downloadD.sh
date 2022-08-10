
#Download  GSM551900 data (LPS stimulation of macrophages)
conda activate pysradb
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332559
pysradb download -p SRP332559 --out-dir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332559/
pysradb metadata SRP332559 > /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332559/sample_sheet.txt

for file in `ls /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332559/SRP332559/*/*.sra`; do
numLines=$(fastq-dump -X 1 -Z --split-spot $file | wc -l)
#echo $file
if [ $numLines -eq 4 ]
then
    echo "$file is single-end"/
    fastq-dump ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332559/SRP332559/single_end/
else
    echo "$file is paired-end"
    fastq-dump --split-files ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332559/SRP332559/paired_end/
fi
done
