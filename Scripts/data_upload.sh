#create folders for data upload
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/data_upload
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/data_upload/rnaseq

cd /omics/groups/OE0219/internal/tinat/raw_data_repo/data_upload/rnaseq
cp /omics/odcf/project/OE0219/tinat/sequencing/rna_sequencing/view-by-pid/OE0219_TINAT_H1299_H2/*/paired/run*/sequence/*fastq.gz ./
cp /omics/odcf/project/OE0219/tinat/sequencing/rna_sequencing/core/run{210622,220202,220217}*/*fastq.gz ./
cp /omics/groups/OE0219/internal/tinat/raw_data_repo/nanopore/EMBL_nanopore_DACandSB939.fastq ./


#copy processed data
cp c010-datasets/Internal/Joschka/data\ upload/RNAseq/*  /omics/groups/OE0219/internal/tinat/raw_data_repo/data_upload/rnaseq/
cp c010-datasets/Internal/Joschka/data\ upload/Riboseq/*  /omics/groups/OE0219/internal/tinat/raw_data_repo/data_upload/riboseq/

#upload to geo
lftp
set ftp:use-hftp false
open ftp://geoftp:rebUzyi1@ftp-private.ncbi.nlm.nih.gov
cd uploads/joschkahey_M7y9GAcU
mirror -R /omics/groups/OE0219/internal/tinat/raw_data_repo/data_upload/rnaseq/
#--> notify geo on webpage +spezifiziere data_upload





