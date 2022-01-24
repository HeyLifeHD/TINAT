#conda activate R4.0
#R
library("megadepth")

##get bigwig file paths
bw_paths <- list.files("/omics/groups/OE0219/internal/tinat/GTEX/data/", pattern=".bw",full.names=TRUE)

## Next, we locate the example annotation BED file
annotation <- "/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.gtf"

## Now we can compute the coverage
bw_cov <- get_coverage(example_bw, op = "mean", annotation = annotation)
bw_cov