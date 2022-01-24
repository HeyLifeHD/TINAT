#Download  GSE140610 data (LPS stimulation of macrophages)
conda activate pysradb
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332555
pysradb download -p SRP332555 --out-dir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332555/
pysradb metadata SRP332555 > /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/v/sample_sheet.txt

for file in `ls /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332555/y/*/*.sra`; do
numLines=$(fastq-dump -X 1 -Z --split-spot $file | wc -l)
#echo $file
if [ $numLines -eq 4 ]
then
    echo "$file is single-end"/
    fastq-dump ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332555/SRP332555/single_end/
else
    echo "$file is paired-end"
    fastq-dump --split-files ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP332555/SRP332555/paired_end/
fi
done

