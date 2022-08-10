#install telescope
conda create --name telescope
conda activate telescope
conda install -c bioconda telescope

#install rediscoverTE

cd /home/heyj/toools
wget http://research-pub.gene.com/REdiscoverTEpaper/software/REdiscoverTE_1.0.1.tar.gz
tar -xvf REdiscoverTE_1.0.1.tar.gz 
rm REdiscoverTE_1.0.1.tar.gz 

conda activate salmon
conda install -c conda-forge r-base=3.4.3 r-essentials
