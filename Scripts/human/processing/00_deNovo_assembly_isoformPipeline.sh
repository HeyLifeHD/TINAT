#run de novo assembly using isoform pipeline
conda activate snakemake

cd  tools/pipeline-nanopore-ref-isoforms

snakemake --config /omics/groups/OE0219/internal/tinat/210709_nanopore_ashish_stringtie_assembly/config.yml

#create bigwig tracks
conda activate deeptools
mkdir /omics/groups/OE0219/internal/tinat/210709_nanopore_ashish_stringtie_assembly/pipeline-nanopore-ref-isoforms/bigwig/
bamCoverage -b /omics/groups/OE0219/internal/tinat/210709_nanopore_ashish_stringtie_assembly/pipeline-nanopore-ref-isoforms/alignments/reads_aln_sorted.bam \
        -o /omics/groups/OE0219/internal/tinat/210709_nanopore_ashish_stringtie_assembly/pipeline-nanopore-ref-isoforms/bigwig/reads_aln_sorted.bw \
        -p 15