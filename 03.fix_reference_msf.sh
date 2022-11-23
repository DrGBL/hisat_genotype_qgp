#path of the imgt-hla alignments
pathAlign=/project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/hisatgenotype_db_original/HLA/alignments/
#path to the orignal fasta files from imgt-hla
pathFasta=/project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/hisatgenotype_db_original/HLA/fasta/
#working directory (will create folders there)
pathWork=/project/richards/guillaume.butler-laporte/bin/hisat-folder/testing_fasta/

cd ${pathWork}

mkdir -p munged
mkdir -p separated
mkdir -p nuc_separated
mkdir -p final

#for 3.50.0
ls ${pathAlign} | grep nuc | sed 's|_nuc.txt||g' | awk '!/DRB/ && !/ClassI/' > list_genes.txt

#for test
#ls ${pathAlign} | grep nuc | sed 's|_nuc.txt||g' | awk '!/DRB/ && !/ClassI/' | awk '/DMA/'> list_genes.txt

#msf file types to build
file_types=("_nuc" "_gen")
#file_types=("_gen")

#loop through that list of genes.
while read line
do
  echo $line
  
  #separate the gen alignment file by chunks of sequences, and remove the header and tail of each alignment. So you're left with only sequences.
  #also replace "*" by "." ("*" means that the sequence is unknown here, but for msf you can replace by ".".
  head -n -3 ${pathAlign}${line}_gen.txt | \
    tail -n +10 | \
    awk '!/gDNA/' | 
    awk '!/cDNA/' | 
    awk '!/AA codon/' | \
    sed 's/^[[:space:]]*|//g' | \
    sed 's|\*|.|g' | \
    sed "s|${line}\.0|${line}\*0|g" | \
    sed "s|${line}\.1|${line}\*1|g" | \
    sed "s|${line}\.2|${line}\*2|g" | \
    sed "s|${line}\.3|${line}\*3|g" | \
    sed "s|${line}\.4|${line}\*4|g" | \
    sed "s|${line}\.5|${line}\*5|g" | \
    sed "s|${line}\.6|${line}\*6|g" | \
    sed "s|${line}\.7|${line}\*7|g" | \
    sed "s|${line}\.8|${line}\*8|g" | \
    sed "s|${line}\.9|${line}\*9|g" | \
    sed 's|^[[:space:]]*$||g' | \
    awk -v RS= -v line=$line '{print > ("separated/" line "_gen_part_" NR ".txt")}' -
    #remove the empty rows
  sed -i '/^[[:space:]]*$/d' separated/${line}_gen*    
    
  #number of files we now have for this alignment, to loop over.
  num_file=$(ls separated | grep ^${line}_gen | wc -l)

  #loop over every alignment subfile, replacing "-" with the reference nucleotide
  for ((file=1; file<=$num_file; file++))
  do
    char_to_replace=$(awk '/^[ 0-9A-Z]*\*[0-9: ]*/' separated/${line}_gen_part_1.txt | head -n 1 | sed -E 's|(^[ 0-9A-Z]*\*[0-9: ]*)[ACTG| \.]*|\1|g' | wc -c)
    #loop over each character in the sequence (sequence starts at character 20, and ends much later (160 is not precise, I just needed a big number).
    #for pos in {20..150}
    for ((pos=$char_to_replace; pos<=160; pos++))
    do
      #extract the reference base pair at that position
      char=$(head -n 1 separated/${line}_gen_part_${file}.txt | head -c ${pos} | tail -c 1)
        
      #now for each line, at this character, if you see "-", replace by the reference
      awk -v pos="$pos" -v char="$char" 'BEGIN{FS=OFS=""} {
        if (substr($0, pos, 1)~"-")
        print(substr($0, 1, pos-1), char, substr($0, pos+1, 150));
        else
        print $0;
      }' separated/${line}_gen_part_${file}.txt > separated/${line}_gen_part_${file}_tmp.txt
      #needed an intermediate file, now I replace it
      mv separated/${line}_gen_part_${file}_tmp.txt separated/${line}_gen_part_${file}.txt
    done
  done


  #now paste 
    
  cat separated/${line}_gen_part_1.txt > munged/${line}_gen_tmp.msf
  mv munged/${line}_gen_tmp.msf munged/${line}_gen.msf
   
  #this big block is due to the weird alleles that are seen at the end of some files. 
  num_rows_one=$(cat separated/${line}_gen_part_1.txt | wc -l )
  for ((file=2; file<=$num_file; file++))
  do
    num_rows_n=$(cat separated/${line}_gen_part_${file}.txt | wc -l)
    if [[ $num_rows_one == $num_rows_n ]]; then
      cut -c 20- separated/${line}_gen_part_${file}.txt | paste munged/${line}_gen.msf - > munged/${line}_gen_tmp.msf
      mv munged/${line}_gen_tmp.msf munged/${line}_gen.msf;
    else
      echo "Last file contains less alleles."
      readarray -t alleles_in_small_doc < <(awk '{print $1}' separated/${line}_gen_part_${file}.txt)
      start=0
      allele_num=0
      for allele in "${alleles_in_small_doc[@]}"
      do
        #echo "test1"
        allele_num=$((allele_num+1))
          
        row_num=$(awk -v allele=$allele '$1==allele {print NR}' separated/${line}_gen_part_1.txt)
        #echo "row_num"
        #echo $row_num
        #echo "start"
        #echo $start
        head -n $((row_num-1)) separated/${line}_gen_part_1.txt | awk -v start=$((start)) 'NR > start { print }' | sed 's|[ACTG\*]|.|g' >> munged/${line}_gen_tmp.msf
        #echo "test2"
        head -n $allele_num separated/${line}_gen_part_${file}.txt | tail -n 1 >> munged/${line}_gen_tmp.msf
          
        start=$row_num
          
        #echo "munged/${line}_gen_tmp.msf"
        #cat munged/${line}_gen_tmp.msf
          
        if [[ $allele_num ==  $num_rows_n && $num_rows_one -gt $start ]]; then
          #echo "dummy";
          awk -v start=$((row_num)) 'NR > start { print }' separated/${line}_gen_part_1.txt | sed 's|[ACTG\*]|.|g' >> munged/${line}_gen_tmp.msf
          #cat munged/${line}_gen_tmp.msf
        fi;
          
      done
        
      #cat munged/${line}_gen_tmp.msf > munged/${line}_gen_tmp_to_play.msf
        
      sed -i 's|[[:space:]]*$||g' munged/${line}_gen_tmp.msf

      cut -c 20- munged/${line}_gen_tmp.msf | paste munged/${line}_gen.msf - > munged/${line}_gen_tmp2.msf
      mv munged/${line}_gen_tmp2.msf munged/${line}_gen.msf
      rm munged/${line}_gen_tmp.msf;
        
      #now pad the result or trim result
        
      n_char_1=$(head -n 1 separated/${line}_gen_part_1.txt | wc -c)
      n_char_n=$(head -n 1 separated/${line}_gen_part_${file}.txt | wc -c)
        
      #trim
      if [[ n_char_1 -gt n_char_n ]]; then
        #echo "tmp"
        sed -E 's|([ACGT\.\*])\s+([ACGT\.\*])|\1\2|g' munged/${line}_gen.msf | sed -E 's|([ACGT\.])\s+([ACGT\.])|\1\2|g' > munged/${line}_gen_tmp.msf
        mv munged/${line}_gen_tmp.msf munged/${line}_gen.msf
        keep_char=$(awk 'NR==1 || length<len {len=length; line=$0} END {print len}' munged/${line}_gen.msf)
        cat munged/${line}_gen.msf > munged/${line}_gen_tmp.msf
        awk -v keep_char=$keep_char '{print(substr($0, 1, keep_char))}' munged/${line}_gen_tmp.msf > munged/${line}_gen.msf
        #mv munged/${line}_gen_tmp2.msf munged/${line}_gen.msf
        rm munged/${line}_gen_tmp.msf;
      fi;
        
      #pad
      if [[ n_char_n -gt n_char_1 ]]; then
        #echo "test more in file 1"
        #sed -i 's|[[:space:]]*$||g' munged/${line}_gen.msf
        sed -E 's|([ACGT\.\*])\s+([ACGT\.\*])|\1\2|g' munged/${line}_gen.msf | sed -E 's|([ACGT\.])\s+([ACGT\.])|\1\2|g' > munged/${line}_gen_tmp.msf
        mv munged/${line}_gen_tmp.msf munged/${line}_gen.msf
        cat munged/${line}_gen.msf > munged/${line}_gen_wtf.msf
        min_char=$(awk 'NR==1 || length<len {len=length; line=$0} END {print len}' munged/${line}_gen.msf)
        max_char=$(awk ' { if ( length > x ) { x = length } }END{ print x }' munged/${line}_gen.msf)
        diff_char=$((max_char-min_char))
        char_to_insert=$(printf '.%.0s' $(seq 1 $diff_char))


        awk -v ct=$char_to_insert -v mc=$max_char 'BEGIN{FS=""; OFS=""} {
          if(length($0) < mc)
          print $0, ct
          else
          print $0
        }' munged/${line}_gen.msf > munged/${line}_gen_tmp.msf
          
        mv munged/${line}_gen_tmp.msf munged/${line}_gen.msf
      fi;      
    fi;
  done
    
  sed -E 's|([ACGT\.\*])\s+([ACGT\.\*])|\1\2|g' munged/${line}_gen.msf | sed -E 's|([ACGT\.])\s+([ACGT\.])|\1\2|g' > munged/${line}_gen_tmp.msf
  mv munged/${line}_gen_tmp.msf munged/${line}_gen.msf

  #now separate the file above in introns and exons
  awk -v line=$line 'BEGIN { FS="|" } { for (i=1; i<=NF; i++) {print $i > ("nuc_separated/" line  "_nuc_part_" i ".txt")}}' munged/${line}_gen.msf
  
  #remove the introns from file 1 (keep only the name of the alleles)
  sed -E 's|(^[[:space:]]*[A-Z0-9]*\*[0-9:A-Z]*[[:space:]]*)[ACTG\.]*|\1|g' nuc_separated/${line}_nuc_part_1.txt > munged/${line}_nuc.msf
  
  num_nuc_files=$(ls nuc_separated | grep ^${line}_nuc | wc -l)
  num_loop=$(( (num_nuc_files-1 )  / 2 ))
  
  for ((i=1; i<=$num_loop; i++))
  do
    j=$((2*i))
    
    paste -d'\0' munged/${line}_nuc.msf nuc_separated/${line}_nuc_part_${j}.txt > munged/${line}_nuc_tmp.msf
    mv munged/${line}_nuc_tmp.msf munged/${line}_nuc.msf
  done



  #now properly align the fasta and msf like real fasta and msf files
  for type in "${file_types[@]}"
  do
    cat munged/${line}${type}.msf | tr ' ' '\n' | sed '/^[[:space:]]*$/d' | sed "s|^${line}\*|>${line}\*|g" | sed 's/|//g' | sed 's|:|x|g' | sed 's|\*|y|g' > final/${line}${type}_tmp.fasta
    seqret -sequence final/${line}${type}_tmp.fasta -outseq final/${line}${type}_tmp2.fasta -osformat2 fasta
    seqret -sequence final/${line}${type}_tmp2.fasta -outseq final/${line}${type}.msf -osformat2 msf

    sed 's|\.||g' final/${line}${type}_tmp.fasta > final/${line}${type}_tmp3.fasta
    seqret -sequence final/${line}${type}_tmp3.fasta -outseq final/${line}${type}.fasta -osformat2 fasta

    sed -i 's|x|:|g'  final/${line}${type}.fasta
    sed -i 's|y|\*|g' final/${line}${type}.fasta
    sed -i 's|-|\.|g' final/${line}${type}.fasta
    sed -i 's|x|:|g'  final/${line}${type}.msf
    sed -i 's|y|\*|g' final/${line}${type}.msf
    sed -i 's|~|\.|g' final/${line}${type}.msf
    sed -i "s|final/${line}${type}\.msf||" final/${line}${type}.msf
    sed -i 's|^>|>HLA:HLA99999 |g'  final/${line}${type}.fasta 
    awk '! /^[0-9 ]+$/' final/${line}${type}.msf > final/${line}${type}_tmp.msf
    mv final/${line}${type}_tmp.msf final/${line}${type}.msf
    rm final/${line}${type}_tmp*
    
    #one little change to the fasta files
    sed 's| |__|g' final/${line}${type}.fasta > final/${line}${type}_tmp.fasta
    grep ">" ${pathFasta}${line}_gen.fasta | sed 's|>||g' | sed 's| |__|g' > final/${line}_list_true_alleles_gen.txt
    grep ">" final/${line}${type}.fasta | sed 's|>||g' | sed 's| |__|g'  > final/${line}_list_wrong_alleles${type}.txt
    paste final/${line}_list_wrong_alleles${type}.txt final/${line}_list_true_alleles_gen.txt > final/${line}_list_alleles${type}.txt
    awk '
    FNR==NR{
      a[$1]=$2
      next
    }
    ($2 in a) && /^>/{
      print ">"a[$2]
      next
    }
    1
    ' final/${line}_list_alleles${type}.txt FS="[> ]"  final/${line}${type}_tmp.fasta > final/${line}${type}.fasta
    sed -i 's|__| |g' final/${line}${type}.fasta
    rm final/${line}${type}_tmp.fasta;
  done;


  #and you're done

done < list_genes.txt

#optional, only if you know where the files need to go already
#cp final/*.fasta /project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/hisatgenotype_db/HLA/fasta/
#cp final/*.msf /project/richards/guillaume.butler-laporte/bin/hisat-folder/hisat_index/hisatgenotype_db/HLA/msf/
