#tested on a python 3.7 conda environment, since they wrote it for 3.7. But probably works for 3.8 or 3.9.
#also require git-lfs (https://git-lfs.github.com/ or with conda install git-lfs)
#and finally need some fastq files, it doesn't matter where they're from, as they're only used to "trick" hisat into building the reference. Here I call them R_1.fastq and R_2.fastq
conda activate hisat-env

cd /some/directory/of/your/choice/

#install hisat-genotype and hisat2
git clone --recurse-submodules https://github.com/DaehwanKimLab/hisat-genotype ./hisatgenotype
echo '{"sanity_check": false}' > hisatgenotype/devel/settings.json
cd hisatgenotype/hisat2
make

cd ../..

#change paths to have the hisat-genotype and hisat2 folders in your path
export PATH=$PWD/hisatgenotype:$PWD/hisatgenotype/hisat2:$PATH
export PYTHONPATH=$PWD/hisatgenotype/hisatgenotype_modules:$PYTHONPATH

#firsrt to test the installation, will test with the pre-compiled reference
#first download it
cd hisatgenotype
mkdir indicies
cd indicies

# genotype_genome
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/genotype_genome_20180128.tar.gz
tar xvzf genotype_genome_20180128.tar.gz
rm genotype_genome_20180128.tar.gz

#grch38
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
tar xvzf grch38.tar.gz
rm grch38.tar.gz
hisat2-inspect grch38/genome > genome.fa
samtools faidx genome.fa

#HISATgenotpye Database
git clone https://github.com/DaehwanKimLab/hisatgenotype_db.git
cd ../..

#now the actual test
hisatgenotype --base hla \
              --threads 5 \
              -z hisatgenotype/indicies \
              --keep-alignment -v --keep-extract \
              -1 R_1.fastq \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat
              
#if the previous code chunk didn't give you HLA calls, then there's something wrong with the installation
              
#now build own reference
#I have managed to make it work for 3.42.0 (below), but it bugs for 3.50.0, 3.49.0, and 3.48.0. Now too sure why.
mkdir -p hisat_index/hisatgenotype_db/
cd hisat_index/hisatgenotype_db/
git clone https://github.com/ANHIG/IMGTHLA HLA
cd HLA
git checkout 3420  # change this to the version you want
git lfs install
git lfs pull
cd ../../..

hisatgenotype --base hla \
              --threads 5 \
              -z hisat_index \
              --keep-alignment -v --keep-extract \
              -1 R_1.fastq \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat
