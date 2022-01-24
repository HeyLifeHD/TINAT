#Download  GSE140610 data (LPS stimulation of macrophages)
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103
for i in {059..100}
  do
    echo SRR3356$i
    fastq-dump --split-3 --gzip SRR3356$i -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103
done

conda activate pysradb
pysradbq metadata SRP073103 > /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/sample_sheet.txt


#merge samples 
mkdir merged_fastq
cd merged_fastq
zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356100.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356099.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356098.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356097.fastq.gz > U-87_DAC_24-4-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356096.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356095.fastq.gzd > U-87_DAC_22-3-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356094.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356093.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356092.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356091.fastq.gz > U-87_24-4-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356090.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356089.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356088.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356087.fastq.gz > U-87_22-3-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356086.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356085.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356084.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356083.fastq.gz > T98G_DAC_29-4-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356082.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356081.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356080.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356079.fastq.gz > T98G_DAC_15-3-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356078.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356077.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356076.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356075.fastq.gz > T98G_29-4-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356074.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356073.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356072.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356071.fastq.gz > T98G_15-3-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356070.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356069.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356068.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356067.fastq.gz > LNT-229_DAC_20-4-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356066.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356065.fastq.gz > LNT-229_DAC_5-2-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356064.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356063.fastq.gz > LNT-229_20-4-15.fastq

zcat /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356062.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356061.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356060.fastq.gz \
/omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRR3356059.fastq.gz > LNT-229_5-2-15.fastq

gzip *

# pysradb download -p SRP073103 --out-dir /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/
# pysradb download -p SRX1689931
# #extract data
# conda deactivate
# mkdir  /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRP073103/single_end
# mkdir  /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRP073103/paired_end

# for file in `ls /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRP073103/*/*.sra`; do
# numLines=$(fastq-dump -X 1 -Z --split-spot $file | wc -l)
# #echo $file
# if [ $numLines -eq 4 ]
# then
#     echo "$file is single-end"/
#     fastq-dump ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRP073103/single_end/
# else
#     echo "$file is paired-end"
#     fastq-dump --split-files ${file} -O /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/SRP073103/SRP073103/paired_end/
# fi
# done

# fastq-dump --split-3 --gzip SRR3356059

