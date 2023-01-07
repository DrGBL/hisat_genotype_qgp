pathHLADat=hisatgenotype_db_original/HLA/hla.dat

mkdir -p hla_dat_file
rm -f hla_dat_file/hla_tmp.dat
#first fix some problems with the hla.dat files.
#note that the hla.dat file is not wrong, it's just that it needs a tweak in order to work with hisat-genotype
#this needs to be done for genes S, T, W. They're first exons are listed as different from 1 in hla.dat for all their alleles (and specifically for the reference allele). This breaks hisat-genotype.


#rm hla_dat_file/hla_tmp.dat

modify="no"
IFS=''

while read line
do

  line_mod=$line
  
  if [[ $line =~ ^DE && $line =~ HLA-P ]]; then
    modify="yes"
    which_gene=P;
  fi
  
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
    
    if [[ $which_gene == P ]]; then
      line_mod=$(echo $line | sed 's|/number=\"3\"|/number=\"1\"|g' | sed 's|/number=\"4\"|/number=\"2\"|g' | sed 's|/number=\"5\"|/number=\"3\"|g' | sed 's|/number=\"6\"|/number=\"4\"|g' | sed 's|/number=\"7\"|/number=\"5\"|g' | sed 's|/number=\"8\"|/number=\"6\"|g' );
    fi
    
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

done < ${pathHLADat}

mv hla_dat_file/hla_tmp.dat hla_dat_file/hla.dat
cp hla_dat_file/hla.dat hisatgenotype_db/HLA/
