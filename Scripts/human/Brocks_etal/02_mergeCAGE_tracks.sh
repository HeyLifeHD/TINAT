#merge CAGE seq data from different strands
conda deactivate
conda activate kentUtils
mkdir /omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE
cd /omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE

#DMSO
bigWigMerge H1299_DMSO_minusS.bw H1299_DMSO_plusS.bw DMSO.bedGraph
bedSort DMSO.bedGraph DMSO.bedGraph
awk 'OFS="\t" {$1="chr"$1; print}' DMSO.bedGraph > DMSO_chr.bedGraph
awk '$1 !~ /GL/' DMSO_chr.bedGraph > DMSO_chr_sub.bedGraph
awk '$1 !~ /MT/' DMSO_chr_sub.bedGraph > DMSO.bedGraph
bedGraphToBigWig DMSO.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DMSO.bigWig

#SB939
bigWigMerge H1299_SB939_minusS.bw H1299_SB939_plusS.bw SB939.bedGraph
bedSort SB939.bedGraph SB939.bedGraph
awk 'OFS="\t" {$1="chr"$1; print}' SB939.bedGraph > SB939_chr.bedGraph
awk '$1 !~ /GL/' SB939_chr.bedGraph > SB939_chr_sub.bedGraph
awk '$1 !~ /MT/' SB939_chr_sub.bedGraph > SB939.bedGraph
bedGraphToBigWig SB939.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     SB939.bigWig

#DAC
bigWigMerge H1299_DAC_minusS.bw H1299_DAC_plusS.bw DAC.bedGraph
bedSort DAC.bedGraph DAC.bedGraph
awk 'OFS="\t" {$1="chr"$1; print}' DAC.bedGraph > DAC_chr.bedGraph
awk '$1 !~ /GL/' DAC_chr.bedGraph > DAC_chr_sub.bedGraph
awk '$1 !~ /MT/' DAC_chr_sub.bedGraph > DAC.bedGraph
bedGraphToBigWig DAC.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DAC.bigWig

#DACSB
bigWigMerge H1299_DACSB_minusS.bw H1299_DACSB_plusS.bw DACSB.bedGraph
bedSort DACSB.bedGraph DACSB.bedGraph
awk 'OFS="\t" {$1="chr"$1; print}' DACSB.bedGraph > DACSB_chr.bedGraph
awk '$1 !~ /GL/' DACSB_chr.bedGraph > DACSB_chr_sub.bedGraph
awk '$1 !~ /MT/' DACSB_chr_sub.bedGraph > DACSB.bedGraph
bedGraphToBigWig DACSB.bedGraph \
 /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
     DACSB.bigWig