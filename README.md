# foldseek 
Software suite for searching and clustering protein structures.
Foldseek is a collaboration between the Söding and Steinegger Lab.

<p align="center"><img src="https://github.com/steineggerlab/foldseek/blob/master/.github/foldseek.png" height="250"/></p>

## Version release
Alpha release: July 24, 2021

## Installation

`foldseek` can be used by compiling from source (see below) or downloading a statically compiled version. It requires a 64-bit system. We recommend using a system with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux).

    # pull docker container
    docker pull steineggerlab/foldseek
    # static Linux AVX2 build
    wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static Linux SSE4.1 build
    wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
    wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

Precompiled binaries for other architectures (ARM64, PPC64LE) and very old AMD/Intel CPUs (SSE2 only) are available at [https://mmseqs.com/foldseek](https://mmseqs.com/foldseek).

### Quick start
    
`easy-search` can search single or multiple queries formatted in pdb/mcif format (flat or gz) against a target database (`example/`) of protein structures. It outputs a tab separated file of the alignments (`.m8`) the fields are `query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits`.

    foldseek easy-search example/d1asha_ example/ aln.m8 tmpFolder
    
The output can be customized with the `--format-output` option e.g. `--format-output "query,target,qaln,taln"` returns the query and target accession and the pairwise alignments in tab separated format. 
You can choose many different [output columns](https://github.com/soedinglab/mmseqs2/wiki#custom-alignment-format-with-convertalis).    

The target database can be pre-processed by `createdb`. This make sense if searched multiple times.
 
    foldseek createdb example/ targetDB
    foldseek easy-search example/d1asha_ targetDB aln.m8 tmpFolder
    
Setup the PDB or AlphaFold using the `databases` module.
    
    # pdb  
    foldseek databases PDB pdb tmp 
    # alphafold db
    foldseek databases AlphafoldDb afdb tmp 

    
### Important parameters

    -s                       adjusyesornot the sensitivity to speed trade-off (default: 7.5, high sensitivity: 9.0)
    --alignment-type         0: 3Di Gotoh-Smith-Waterman (fast), 1: TMalign
    -c                       list matches above this fraction of aligned (covered) residues (see --cov-mode) (default: 0.0) 
    --cov-mode               0: coverage of query and target, 1: coverage of target, 2: coverage of query

### Main Modules

* `easy-search`       fast protein structure search  
* `createdb`          create a database from protein structures (PDB,mmCIF, mmJSON)
* `databases`         download pre-assembled databases

### Compile from source

Compiling `foldseek` from source has the advantage of system-specific optimizations, which should improve its performance. To compile it `git`, `g++` (4.9 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the foldseek binary will be located in the `build/bin` directory.

    git clone https://github.com/steineggerlab/foldseek.git
    cd foldseek
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make -j
    make install
    export PATH=$(pwd)/foldseek/bin/:$PATH

:exclamation: If you want to compile `foldseek` on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and `foldseek` will not be able to run multi-threaded. Adjust the `cmake` call above to:

    CC="$(brew --prefix)/bin/gcc-10" CXX="$(brew --prefix)/bin/g++-10" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..


## Hardware requirements

`foldseek` will scale its memory consumption based on the available main memory of the machine. `foldseek` needs a 64-bit CPU with at least the SSE2 instruction set to run.
