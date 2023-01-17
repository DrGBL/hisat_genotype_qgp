#tested on a python 3.7 conda environment, since they wrote it for 3.7. But probably works for 3.8 or 3.9.
#also require git-lfs (https://git-lfs.github.com/ or with conda install git-lfs)
#and finally need some fastq files, it doesn't matter where they're from, as they're only used to "trick" hisat into building the reference. Here I call them R_1.fastq and R_2.fastq
conda activate hisat-env

cd /some/directory/of/your/choice/

export PYTHONUNBUFFERED=true

#install hisat-genotype and hisat2
git clone --recurse-submodules https://github.com/DrGBL/hisat-genotype ./hisatgenotype
echo '{"sanity_check": false}' > hisatgenotype/devel/settings.json
cd hisatgenotype/hisat2
make

cd ../..

#change paths to have the hisat-genotype and hisat2 folders in your path
export PATH=$PWD/hisatgenotype:$PWD/hisatgenotype/hisat2:$PATH
export PYTHONPATH=$PWD/hisatgenotype/hisatgenotypeules:$PYTHONPATH

#download hla database, note that I used 3.50.0 for all genes except for DPB1 due to a strange misalignment in some genes (e.g. DPB1*04:01:01:34)
mkdir -p hisat_index/
cd hisat_index
mkdir -p hisatgenotype_db/
mkdir -p hisatgenotype_db/HLA
mkdir -p hisatgenotype_db/HLA/msf
mkdir -p hisatgenotype_db/HLA/fasta
mkdir -p hisatgenotype_db_3490/
mkdir -p hisatgenotype_db_original/HLA/alignments/
mkdir -p hisatgenotype_db_original/HLA/fasta/

cd hisatgenotype_db_3490/
git clone --single-branch --branch 3490 https://github.com/ANHIG/IMGTHLA HLA
cd ..

#move files around
mv hisatgenotype_db_3490/HLA/hla.dat hisatgenotype_db_original/HLA/
mv hisatgenotype_db_3490/HLA/alignments/* hisatgenotype_db_original/HLA/alignments/
mv hisatgenotype_db_3490/HLA/fasta/* hisatgenotype_db_original/HLA/fasta/

rm -r hisatgenotype_db_3490/


#adjust imgt hla files
sh ../01.fix_hla_dat_file.sh
rm -r hla_dat_file
bash ../02.fix_reference_msf.sh
rm -r final munged nuc_separated separated hisatgenotype_db_original

#download cyp and codis folders from here https://github.com/DaehwanKimLab/hisatgenotype_db
#but note that I haven't updated them at all
svn checkout https://github.com/DaehwanKimLab/hisatgenotype_db/trunk/CODIS hisatgenotype_db/CODIS
svn checkout https://github.com/DaehwanKimLab/hisatgenotype_db/trunk/CYP hisatgenotype_db/CYP

#obtain snps and haplotypes for graph genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz
Rscript ../03.fix_repeat_rsid.R
#hisat2_extract_snps_haplotypes_UCSC.py Homo_sapiens_assembly38_original.fasta snp151Common_fixed.txt.gz snp151Common_fixed


#download hrch38 with the HLA haplotypes
wget https://media.githubusercontent.com/media/broadinstitute/gatk/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz
gzip -d Homo_sapiens_assembly38.fasta.gz

#modify references so that there's only chromosomes and alternative scaffolds (no decoys and other nonsense)
#each chr6 scaffold separately, unfortunately
mv Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38_original.fasta
#chr6
mkdir -p chr6
grep "^>" Homo_sapiens_assembly38_original.fasta | \
  awk '/^>chr[0-9XYM]+/' | \
  awk '!/unlocalized/' | \
  sed -E 's|^>([0-9a-zA-Z_]*).*|\1|g' | \
  awk '/^chr[0-9MXY]*$/' > chr6/list_contigs_to_keep_chr6.txt
seqtk subseq Homo_sapiens_assembly38_original.fasta chr6/list_contigs_to_keep_chr6.txt > chr6/Homo_sapiens_assembly38_chr6.fasta
mkdir -p chr6/grch38
hisat2-build -p 30 chr6/Homo_sapiens_assembly38_chr6.fasta chr6/grch38/genome #-p is number of threads
hisat2-inspect chr6/grch38/genome > chr6/genome.fa
samtools faidx chr6/genome.fa
rm chr6/list_contigs_to_keep_chr6.txt chr6/Homo_sapiens_assembly38_chr6.fasta

#chr6_GL00025*
for n in chr6_GL000250 chr6_GL000251 chr6_GL000252 chr6_GL000253 chr6_GL000254 chr6_GL000255 chr6_GL000256
do
  mkdir -p ${n}
  grep "^>" Homo_sapiens_assembly38_original.fasta | \
    awk '/^>chr[0-9XYM]+/' | \
    awk '!/unlocalized/' | \
    sed -E 's|^>([0-9a-zA-Z_]*).*|\1|g' | \
    awk -v n=$n '/^chr[0-9MXY]*$/ || $0 ~ n' | \
    awk '!/^chr6$/' > ${n}/list_contigs_to_keep_${n}.txt
  seqtk subseq Homo_sapiens_assembly38_original.fasta ${n}/list_contigs_to_keep_${n}.txt > ${n}/Homo_sapiens_assembly38_${n}.fasta
  mkdir -p ${n}/grch38
  hisat2-build -p 30 ${n}/Homo_sapiens_assembly38_${n}.fasta ${n}/grch38/genome #-p is number of threads
  hisat2-inspect ${n}/grch38/genome > ${n}/genome.fa
  samtools faidx ${n}/genome.fa
  rm ${n}/list_contigs_to_keep_${n}.txt ${n}/Homo_sapiens_assembly38_${n}.fasta
done

#build genome_genotype
for n in chr6 chr6_GL000250 chr6_GL000251 chr6_GL000252 chr6_GL000253 chr6_GL000254 chr6_GL000255 chr6_GL000256
do
  cp -r hisatgenotype_db ${n}/
  cd ${n}
  hisat2_extract_snps_haplotypes_UCSC.py genome.fa ../snp151Common_fixed.txt.gz snp151Common_fixed_${n}
  python ../../hisatgenotype/hisatgenotype_tools/hisatgenotype_build_genome_mod.py --base genotype_genome --locus-list codis,cyp,hla --prefix_snps snp151Common_fixed_${n} > hisat_index_${n}.log 2>&1
  hisat2-build -p 30 --snp hla.index.snp --haplotype hla.haplotype hla_backbone.fa hla.graph > hisat_hla_index_${n}.log 2>&1
  rm snp151Common_fixed_${n}*
  cd ..
done

#build genotype_genome graph genome
#the following would be the code, but I did this on dna-nexus since it requires >200GB of ram
#hisat2-build -p 96 --snp genotype_genome.index.snp --haplotype genotype_genome.haplotype genotype_genome.fa genotype_genome

#hence, put the genotype_genome files in tar balls for export to dnanexus
for n in chr6 chr6_GL000250 chr6_GL000251 chr6_GL000252 chr6_GL000253 chr6_GL000254 chr6_GL000255 chr6_GL000256
do
  cd ${n}
  tar -czvf genotype_genome_${n}.tar.gz genotype_genome*
  cd ..
done

#now run hisat-genotype!
cd ..

hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6 \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list A,B,C,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DPB2,DQA1,DQA2,DQB1,DRA,DRB1,DRB5,E,F,G,H,HFE,J,K,L,MICA,MICB,N,P,S,T,TAP1,TAP2,U,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6 > hisat_build_chr6.log 2>&1

#chr6_GL000250
hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6_GL000250/ \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list A,DMA,DMB,DOB,DPA1,DPA2,DPB1,DPB2,DQA2,DRB3,G,H,J,K,MICB,N,P,T,TAP1,TAP2,U,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6_GL000250/ > hisat_build_chr6_GL000250.log 2>&1

#chr6_GL000251
hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6_GL000251 \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list A,B,C,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DPB2,DQA1,DQA2,DQB1,DRA,DRB1,DRB3,E,F,G,H,J,K,L,MICA,MICB,N,P,S,T,TAP1,TAP2,U,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6_GL000251 > hisat_build_chr6_GL000251.log 2>&1

#chr6_GL000252
hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6_GL000252 \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list A,C,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DQA1,DQA2,DQB1,DRA,DRB1,E,F,G,H,J,K,L,MICB,N,P,T,TAP1,TAP2,U,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6_GL000252 > hisat_build_chr6_GL000252.log 2>&1

#chr6_GL000253
hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6_GL000253 \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list A,B,C,DMA,DMB,DOA,DPA1,DPA2,DPB1,DPB2,DQA1,DQA2,DQB1,DRA,DRB1,DRB4,E,F,G,H,J,K,L,MICA,N,P,S,T,TAP1,TAP2,U,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6_GL000253 > hisat_build_chr6_GL000253.log 2>&1

#chr6_GL000254
hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6_GL000254 \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list B,C,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DQA1,DQA2,DQB1,DRA,DRB4,E,G,H,J,MICB,P,T,TAP2,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6_GL000254 > hisat_build_chr6_GL000254.log 2>&1

#chr6_GL000255
hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6_GL000255 \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list A,B,C,DMA,DMB,DOA,DPA1,DPA2,DPB1,DPB2,DQA1,DQA2,DQB1,DRA,DRB1,DRB3,E,F,G,H,J,K,L,MICA,MICB,N,P,S,T,TAP1,TAP2,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6_GL000255 > hisat_build_chr6_GL000255.log 2>&1

#chr6_GL000256
hisatgenotype --base hla \
              --threads 10 \
              -z hisat_index/chr6_GL000256 \
              -v --keep-extract \
              -1 R_1.fastq \
              --locus-list A,B,C,DMA,DMB,DOA,DOB,DPA1,DPA2,DPB1,DPB2,DQA1,DQA2,DQB1,DRA,DRB4,E,F,G,H,J,K,L,MICA,MICB,N,P,S,T,TAP1,TAP2,U,V,W \
              -2 R_2.fastq \
              --out-dir ./tmp_hisat/chr6_GL000256 > hisat_build_chr6_GL000256.log 2>&1









