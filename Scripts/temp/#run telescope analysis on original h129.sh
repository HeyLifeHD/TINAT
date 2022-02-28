#run telescope analysis on original h1299 data
conda activate telescope

mkdir /omics/groups/OE0219/internal/tinat/211007_shortRead_telescope_quantification/
cd /omics/groups/OE0219/internal/tinat/211007_shortRead_telescope_quantification/

telescope assign /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/HISAT2/aligned_sorted/AS-641773-LR-57123.sorted.bam \
/omics/groups/OE0219/internal/repguide/data/repeats/repeats_hg19.gtf \
--ncpu 5 --outdir /omics/groups/OE0219/internal/tinat/211007_shortRead_telescope_quantification/

telescope assign /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/HISAT2/aligned_sorted/AS-641773-LR-57123.sorted.bam \
./transcripts.gtf \
--ncpu 5 --outdir /omics/groups/OE0219/internal/tinat/211007_shortRead_telescope_quantification/