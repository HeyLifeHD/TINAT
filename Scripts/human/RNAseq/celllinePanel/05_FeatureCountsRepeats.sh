##Count overlap or reads with repeats
conda activate featureCounts

files=`ls /omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification/HISAT2/aligned_sorted/*.bam`
#not allowing multimapping
mkdir /omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_FeatureCountsRepeats
featureCounts -p -T 7 -f -F 'GTF' -a  /omics/groups/OE0219/internal/repguide/data/repeats/repeats_hg19.gtf.gz \
-o /omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_FeatureCountsRepeats/repeats_featureCounts.txt $files 
