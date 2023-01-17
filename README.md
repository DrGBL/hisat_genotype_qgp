# hisat_genotype_qgp
This git provides code to create your own indexes for use with HISAT-genotype. It assumes that you are using the modified version of hisat-genotype, available here: https://github.com/DrGBL/hisat-genotype

Note that starting from IMGT-HLA version 3.42.0, HISAT-genotype cannot work with alleles for which the full sequence length is not available. Specifically, the alleles that are present in the _nuc files, but not in the _gen files, cannot be used.

To run the code here, simply put files 01 to 04 in some folder of your choice, adjust the path to that folder at the top of file 04, and run file 04.

This will give 8 folders necessary to run hisat-genotype on all genes appropriately. The way to run it is shown at the bottom of file 04.
