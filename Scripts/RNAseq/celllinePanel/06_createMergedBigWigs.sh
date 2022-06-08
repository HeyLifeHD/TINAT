cd c010-datasets/Internal/Joschka/cell\ line\ pannel\ bw\ renamed/

#rename files
for f in *; do mv "$f" "${f// /_}"; done 
for f in *; do mv "$f" "${f//_+_/_}"; done 

#get all different 
conda activate kentUtils
conditions=`ls -1 ./ | awk -F'_R' '{print $1}' | awk '!a[$0]++'`
for file in echo $conditions;do
    ls ./${file}*
    bigWigMerge ${file}*  ${file}.bedgraph 
    bedSort  ${file}.bedgraph ${file}.bedgraph 
    bedGraphToBigWig ${file}.bedgraph \
        /omics/groups/OE0219/internal/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes \
        ${file}.bigwig 
done