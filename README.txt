Minimac3 is a lower memory and more computationally efficient
implementation of the genotype imputation algorithms in 
minimac and minimac2.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac4 FOR DOCUMENTATION >>>

## CLONE MINIMAC4 AND COMPILE
 
git clone https://github.com/Santy-8128/Minimac4
cd Minimac3/
make

A typical Minimac4 command line for imputation is as follows

../bin/Minimac4 --refHaps refPanel.m3vcf.gz \ 
                --haps targetStudy.vcf.gz \
                --prefix testRun

Here refPanel.m3vcf.gz is the reference panel used in M3VCF format (e.g. 1000 Genomes), 
targetStudy.vcf.gz is the phased GWAS data in VCF format, 
and testRun is the prefix for the output files.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac4 FOR DOCUMENTATION>>>
