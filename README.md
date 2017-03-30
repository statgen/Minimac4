# Minimac4

THIS IS FOR TESTING PURPOSES ONLY !!!

Minimac4 is a lower memory and more computationally efficient
implementation of the genotype imputation algorithms in 
minimac/mininac2/minimac3.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac4 FOR DOCUMENTATION >>>

Users should follow the following steps to compile Minimac4 
(if they downloaded the source files) or should skip them
(if they downloaded the binary executable).

## EXTRACT MINIMAC3 AND COMPILE
 
tar -xzvf Minimac4.v1.tar.gz
cd Minimac4/
make

A typical Minimac4 command line for imputation is as follows

../bin/Minimac4 --refHaps refPanel.m3vcf \ 
                --haps targetStudy.vcf \
                --prefix testRun

Here refPanel.vcf is the reference panel used in M3VCF format (e.g. 1000 Genomes), 
targetStudy.vcf is the phased GWAS data in VCF format, 
and testRun is the prefix for the output files.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac4 FOR DOCUMENTATION>>>
