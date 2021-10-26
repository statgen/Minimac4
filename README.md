# Minimac4

Minimac4 is a lower memory and more computationally efficient
implementation of the genotype imputation algorithms in 
minimac/mininac2/minimac3.

<<< SEE http://genome.sph.umich.edu/wiki/Minimac4 FOR DOCUMENTATION >>>

## Prerequisites
Automatic installation of Minimac4 requires [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget) and cmake >= v3.2.

## Installation
The easiest way to install Minimac4 and its dependencies is to use cget:
```bash
cget install --prefix <install_prefix> statgen/Minimac4
```

Alternatively, you can install manually:
```bash
cd Minimac4
cget install -f ./requirements.txt                      # Install dependencies locally.
mkdir build && cd build                                 # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                    # Build.
make install                                            # Install
```

## Usage
See `minimac4 --help` for detailed usage.

A typical Minimac4 command line for imputation is as follows
```bash
minimac4 reference.msav target.bcf > imputed.bcf
```

Here reference.msav is a reference panel (e.g. 1000 Genomes) compressed with MVCF encoding, 
target.bcf is an indexed BCF containing phased genotype array data, 
and imputed.bcf is the imputed output.
