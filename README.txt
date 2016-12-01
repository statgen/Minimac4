Minimac3 is a lower memory and more computationally efficient
implementation of the genotype imputation algorithms in 
minimac and minimac2.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac3 FOR DOCUMENTATION >>>

Users should follow the following steps to compile Minimac3 
(if they downloaded the source files) or should skip them
(if they downloaded the binary executable).

## EXTRACT MINIMAC3 AND COMPILE
 
tar -xzvf Minimac3.v1.tar.gz
cd Minimac3/
make

A typical Minimac3 command line for imputation is as follows

../bin/Minimac3 --refHaps refPanel.vcf \ 
                --haps targetStudy.vcf \
                --prefix testRun

Here refPanel.vcf is the reference panel used in VCF format (e.g. 1000 Genomes), 
targetStudy.vcf is the phased GWAS data in VCF format, 
and testRun is the prefix for the output files.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac3 FOR DOCUMENTATION>>>
