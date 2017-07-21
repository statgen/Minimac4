# Minimac4

THIS IS FOR TESTING PURPOSES ONLY !!!

Minimac4 is a lower memory and more computationally efficient
implementation of the genotype imputation algorithms in 
minimac/mininac2/minimac3.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac4 FOR DOCUMENTATION >>>

Users should follow the following steps to compile Minimac4 
(if they downloaded the source files) or should skip them
(if they downloaded the binary executable).

## Installation
The easiest way to install Minimac4 and its dependencies is to use [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget).
```bash
cget install --prefix <install_prefix> statgen/savvy
```

Alternatively, you can setup a dev environment cmake directly.
```bash
cd Minimac4
cget install -f ./requirements.txt                      # Install dependencies locally.
mkdir build && cd build                                 # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                    # Build.
```

## Usage
A typical Minimac4 command line for imputation is as follows
```bash
minimac4 --refHaps refPanel.m3vcf \ 
                --haps targetStudy.vcf \
                --prefix testRun
```
Here refPanel.vcf is the reference panel used in M3VCF format (e.g. 1000 Genomes), 
targetStudy.vcf is the phased GWAS data in VCF format, 
and testRun is the prefix for the output files.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac4 FOR DOCUMENTATION>>>
