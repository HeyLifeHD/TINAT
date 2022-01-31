#get data
cd /omics/groups/OE0219/internal/tinat/Encode/data_new/
cat files.txt | parallel --gnu "wget {}"

wget https://www.encodeproject.org/metadata/?type=Experiment&cart=%2Fcarts%2F8f6f12db-ad91-4c8e-be79-aa5e45200034%2F&files.output_category=raw+data 