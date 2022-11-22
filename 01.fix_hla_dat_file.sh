pathAlign=/project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/hisatgenotype_db_original/HLA/alignments/
pathFasta=/project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/hisatgenotype_db_original/HLA/fasta/
pathWork=/project/richards/guillaume.butler-laporte/bin/hisat-folder/testing_fasta/
pathHLADat=/project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/hisatgenotype_db_original/HLA/hla.dat

cd ${pathWork}

mkdir -p munged
mkdir -p separated
mkdir -p final
mkdir -p ref_alleles_only
mkdir -p hla_dat_file

#first fix some problems with the hla.dat files.
#note that the hla.dat file is not wrong, it's just that it needs a tweak in order to work with hisat-genotype
#this needs to be done for genes S, T, W. They're first exons are listed as different from 1 in hla.dat for all their alleles (and specifically for the reference allele). This breaks hisat-genotype.

#it also needs to be fixed to merge HLA-DRB 3-4-5 together as HLA-DRB345

sed 's|DE   HLA-DRB3\*|DE   HLA-DRB345\*399|g' ${pathHLADat} | \
  sed 's|DE   HLA-DRB4\*|DE   HLA-DRB345\*499|g' | \
  sed 's|DE   HLA-DRB5\*|DE   HLA-DRB345\*599|g' > hla_dat_file/hla.dat

rm hla_dat_file/hla_tmp.dat

modify="no"
IFS=''
cat hla_dat_file/hla.dat |
while read line
do

  line_mod=$line

  if [[ $line =~ ^DE && $line =~ HLA-S ]]; then
    modify="yes"
    which_gene=S;
  fi
  
  if [[ $line =~ ^DE && $line =~ HLA-T ]]; then
    modify="yes"
    which_gene=T;
  fi
  
  if [[ $line =~ ^DE && $line =~ HLA-W ]]; then
    modify="yes"
    which_gene=W;
  fi
  
  if [[ $line =~ ^SQ ]]; then
    modify="no";
  fi
  
  if [[ $line =~ ^FT && $modify == "yes" ]]; then
    
    if [[ $which_gene == S ]]; then
      line_mod=$(echo $line | sed 's|/number=\"6\"|/number=\"1\"|g' | sed 's|/number=\"7\"|/number=\"2\"|g' | sed 's|/number=\"8\"|/number=\"3\"|g' );
    fi
    
    if [[ $which_gene == T ]]; then
      line_mod=$(echo $line | sed 's|/number=\"4\"|/number=\"1\"|g' | sed 's|/number=\"5\"|/number=\"2\"|g' | sed 's|/number=\"6\"|/number=\"3\"|g' | sed 's|/number=\"7\"|/number=\"4\"|g' );
    fi
    
    if [[ $which_gene == W ]]; then
      line_mod=$(echo $line | sed 's|/number=\"3\"|/number=\"1\"|g' | sed 's|/number=\"4\"|/number=\"2\"|g' | sed 's|/number=\"5\"|/number=\"3\"|g' | sed 's|/number=\"6\"|/number=\"4\"|g' | sed 's|/number=\"7\"|/number=\"5\"|g' | sed 's|/number=\"8\"|/number=\"6\"|g');
    fi;
    
  fi
  
  echo $line_mod >> hla_dat_file/hla_tmp.dat;

done

mv hla_dat_file/hla_tmp.dat hla_dat_file/hla.dat

#optional, only if you already know where to send hla.dat
#cp hla_dat_file/hla.dat ../hisat_index/hisatgenotype_db/HLA/
