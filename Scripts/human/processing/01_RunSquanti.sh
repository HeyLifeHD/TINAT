 conda activate SQANTI3.env
cd tools/SQANTI3-4.1/ 
export PYTHONPATH=$PYTHONPATH:/home/heyj/tools/cDNA_Cupcake/sequence/

python sqanti3_qc.py \
/omics/groups/OE0219/internal/tinat/210709_nanopore_ashish_stringtie_assembly/pipeline-nanopore-ref-isoforms/results/annotation/str_merged_unspecifiedRemoved.gff \
/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf \
/omics/groups/OE0219/internal/tinat/raw_data_repo/references/hg19.fa \
--min_ref_len 50 \
--force_id_ignore \
--aligner_choice minimap2 \
--polyA_motif_list data/polyA_motifs/mouse_and_human.polyA_motif.txt \
-t 12 \
-o 210715_Integration \
-d /home/heyj/temp/210715.3_integrated_deNovo_analysis_LongReadGFFadapted_ShortReadDACSB_CAGEConsenus/ \
--saturation \
--report both \
--cage_peak /omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/consensus_clusters_SQUANTI.bed \
--short_reads /omics/groups/OE0219/internal/tinat/raw_data_repo/RNA_seq/DMSO_DAC_DACandSB939/DACandSB939.fofn