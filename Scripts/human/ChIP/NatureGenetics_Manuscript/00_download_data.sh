#Download  GSE140610 data (LPS stimulation of macrophages)
conda activate pysradb
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/fastq
#save srr list
#nano srr_list.txt
cd /omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/fastq
#loop over file and download
conda create --name sra_tools
conda activate sra_tools
conda install -c bioconda sra-tools

prefetch --option-file srr_list.txt 
fasterq-dump ./* -O ./
