#start characterizing de novo assembly
#gtf == /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted.gtf
#gff compare == /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare 

#get stat information of gff compare vs reference
cd /omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/ 
##n. assembled genes:
cat gffCompare.mergedTranscripts.gtf.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l
#68526

##n. novel genes:
cat gffCompare.mergedTranscripts.gtf.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l
#1455

##n. novel transcripts:
cat gffCompare.mergedTranscripts.gtf.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l
#2157

##n. transcripts matching annotation:
cat gffCompare.mergedTranscripts.gtf.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l
#207411