mkdir  /omics/groups/OE0219/internal/tinat/manuscript_suppData

#supplementary data 1
mkdir  /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_1
for file in /omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef_analysis/results/PostDE/*/DEG.csv; do 
parent_dir_name=$(echo $file | tr  "\/" "\t" | cut -f10 )
echo $parent_dir_name
cp $file /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_1/${parent_dir_name}_DEG.csv
done

#supplementary data 1
mkdir  /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_2
cp /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_2/

#supplementary data 3
mkdir  /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_3
for file in /omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/results/PostDE/*/DEG.csv; do 
parent_dir_name=$(echo $file | tr  "\/" "\t" | cut -f10 )
echo $parent_dir_name
cp $file /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_3/${parent_dir_name}_DET.csv
done

#supplementary data 4
mkdir  /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_4
cp /omics/groups/OE0219/internal/tinat/Encode/220128_deNovo_quantification_H1299ref_analysis//results/tables/tpm_mean.csv /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_4/

#supplementary data 5
mkdir /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_5
cp /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/transdecoder_default_topStrand_8aa/longest_orfs_validated_induced_DACandSB939vsDMSO_forJens_novel.fa /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_5/

#supplementary data 6
mkdir /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_6
cp /omics/groups/OE0219/internal/tinat/proteomics/data/Limma_analysis_SBDAC_nascent_proteome.txt /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_6/
cp /omics/groups/OE0219/internal/tinat/proteomics/data/Limma_analysis_TINATs_translatome.txt /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_6/

#supplementary data 7
mkdir /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_7
cp /omics/groups/OE0219/internal/tinat/integration/peptidomics/data/220301_Peptide_Results_new.xlsx  /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_7/

#supplementary data 8
mkdir /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_8
cp /omics/groups/OE0219/internal/tinat/220316_cellline_deNovo_assembly/gffCompare.annotated.sorted.gtf /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_8/

#supplementary data 9
mkdir /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_9
cp /omics/groups/OE0219/internal/tinat/Cellline_panel/220319_cellline_deNovo_assembly_quantification_analysis/results/PostDE/DACSB_vs_DMSO/DEG.csv  /omics/groups/OE0219/internal/tinat/manuscript_suppData/supp_data_9/