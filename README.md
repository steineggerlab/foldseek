# Foldseek 
Foldseek enables fast and sensitive comparisons of large structure sets. It reaches sensitivities similar to state-of-the-art structural aligners while being at least 20,000 times faster.

<p align="center"><img src="https://github.com/steineggerlab/foldseek/blob/master/.github/foldseek.png" height="250"/></p>

## Publications

[van Kempen M, Kim S, Tumescheit C, Mirdita M, Söding J, and Steinegger M. Foldseek:  fast and accurate protein structure search. bioRxiv, doi:10.1101/2022.02.07.479398  (2022)](https://www.biorxiv.org/content/10.1101/2022.02.07.479398)

## Webserver 
Search your protein structures against the [AlphaFoldDB](https://alphafold.ebi.ac.uk/) and [PDB](https://www.rcsb.org/) in seconds using our Foldseek webserver:

[🚀search.foldseek.com](https://search.foldseek.com)

## Installation

`foldseek` can be used by compiling from source (see below) or downloading a statically compiled version. It requires a 64-bit system. We recommend using a system with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux).

    # static Linux AVX2 build
    wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static Linux SSE4.1 build
    wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
    wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # conda installer 
    conda install -c conda-forge -c bioconda foldseek

Precompiled binaries for other architectures (ARM64, PPC64LE) and very old AMD/Intel CPUs (SSE2 only) are available at [https://mmseqs.com/foldseek](https://mmseqs.com/foldseek).

### Quick start
    
`easy-search` can search single or multiple queries formatted in PDB/mmCIF format (flat or `.gz`) against a target database (`example/`) of protein structures. It outputs a tab-separated file of the alignments (`.m8`) the fields are `query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits`.

    foldseek easy-search example/d1asha_ example/ aln.m8 tmpFolder
    
The output can be customized with the `--format-output` option e.g. `--format-output "query,target,qaln,taln"` returns the query and target accession and the pairwise alignments in tab separated format. 
You can choose many different [output columns](https://github.com/soedinglab/mmseqs2/wiki#custom-alignment-format-with-convertalis).    

The target database can be pre-processed by `createdb`. This make sense if searched multiple times.
 
    foldseek createdb example/ targetDB
    foldseek easy-search example/d1asha_ targetDB aln.m8 tmpFolder

### Important parameters
    -s                       adjust the sensitivity to speed trade-off.
                             lower is faster, higher more sensitive (fast: 7.5, highest sensitivity (default): 9.5)
    --max-seqs               adjust the amount of prefilter that are handed to the alignment. 
                             Increasing it can lead to more hits (default: 300)
    --alignment-type         0: 3Di Gotoh-Smith-Waterman (local, not recommended), 
                             1: TMalign (global), 
                             2: 3Di+AA Gotoh-Smith-Waterman (local, default)
    -c                       list matches above this fraction of aligned (covered) residues (see --cov-mode) (default: 0.0) 
    --cov-mode               0: coverage of query and target, 1: coverage of target, 2: coverage of query



### Databases 
Setup the PDB or AlphaFoldDB using the `databases` module.
    
    # pdb  
    foldseek databases PDB pdb tmp 
    # alphafold db
    foldseek databases Alphafold/Proteome afdb tmp 

We currently support the following databases: 
```
  Name                  Type            Taxonomy        Url
- Alphafold/Proteome    Aminoacid            yes        https://alphafold.ebi.ac.uk/
- Alphafold/Swiss-Prot  Aminoacid            yes        https://alphafold.ebi.ac.uk/
- PDB                   Aminoacid            yes        https://www.rcsb.org
```
    

### Main Modules

* `easy-search`       fast protein structure search  
* `createdb`          create a database from protein structures (PDB,mmCIF, mmJSON)
* `databases`         download pre-assembled databases

### TMalign/TMscore 
Foldseek supports to realign hits using TMalign as well as rescoring alignments using TMscore. 

In case of the alignment type (`--alignment-type 1`) tmalign we sort the results by the TMscore normalized by query length. We write the TMscore into the e-value(=TMscore) as well as into the score(=TMscore*100) field.

```
foldseek easy-search example/d1asha_ example/ aln tmp --alignment-type 1
```

It is possible to compute TMscores for the kind of alignment output (e.g. 3Di/AA) using the following commands: 
```
foldseek createdb example/ targetDB
foldseek createdb example/ queryDB
foldseek search queryDB targetDB aln tmpFolder -a
foldseek aln2tmscore queryDB targetDB aln aln_tmscore
foldseek createtsv queryDB targetDB aln_tmscore aln_tmscore.tsv
```

In the output is the query and target identifier, TMscore, translation(3) and rotation vector=(3x3) (`query,target,TMscore,t[0-2],u[0-2][0-2]`)

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

:exclamation: If you want to compile `foldseek` on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP (by default) and `foldseek` will not be able to run multi-threaded. Adjust the `cmake` call above to:

    CC="$(brew --prefix)/bin/gcc-11" CXX="$(brew --prefix)/bin/g++-11" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..


## Hardware requirements

`foldseek` will scale its memory consumption based on the available main memory of the machine. `foldseek` needs a 64-bit CPU with at least the SSE2 instruction set to run.
