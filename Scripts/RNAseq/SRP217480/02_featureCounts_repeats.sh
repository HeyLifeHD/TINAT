##Count overlap or reads with repeats
conda activate featureCounts

files=`ls /omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/HISAT2/aligned_sorted/*.bam`
#not allowing multimapping
featureCounts -p -T 5 -f -F 'GTF' -a /omics/groups/OE0219/internal/repguide/data/repeats/repeats_hg19.gtf.gz \
-o /omics/groups/OE0219/internal/tinat/210928_SRP217480_processing_knownRef/repeats_featureCounts.txt $files 
