# Foldseek 
Foldseek enables fast and sensitive comparisons of large structure sets.

<p align="center"><img src="https://github.com/steineggerlab/foldseek/blob/master/.github/foldseek.png" height="250"/></p>

## Publications
[van Kempen M, Kim S, Tumescheit C, Mirdita M, Lee J, Gilchrist C, Söding J, and Steinegger M. Foldseek: fast and accurate protein structure search. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)](https://www.nature.com/articles/s41587-023-01773-0)

[Barrio-Hernandez I, Yeo J, Jänes J, Wein T, Varadi M, Velankar S, Beltrao P and Steinegger M. Clustering predicted structures at the scale of the known protein universe. biorxiv, doi:10.1101/2023.03.09.531927 (2023)](https://www.biorxiv.org/content/10.1101/2023.03.09.531927v1)
# Table of Contents

- [Foldseek](#foldseek)
- [Webserver](#webserver)
- [Installation](#installation)
- [Memory requirments](#memory-requirments)
- [Tutorial Video](#tutorial-video)
- [Documentation](#documentation)
- [Quick Start](#quick-start)
  - [Search](#search)
    - [Output](#output-search)
    - [Important Parameters](#important-search-parameters)
    - [Alignment Mode](#alignment-mode)
  - [Databases](#databases)
    - [Create Custom Databases and Indexes](#create-custom-databases-and-indexes)
  - [Cluster](#cluster)
    - [Output](#output-cluster)
    - [Important Parameters](#important-cluster-parameters)
- [Main Modules](#main-modules)
- [Examples](#examples)

## Webserver 
Search your protein structures against the [AlphaFoldDB](https://alphafold.ebi.ac.uk/) and [PDB](https://www.rcsb.org/) in seconds using our Foldseek webserver: [search.foldseek.com](https://search.foldseek.com) 🚀

## Installation
```
# Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

# Linux SSE4.1 build (check using: cat /proc/cpuinfo | grep sse4_1)
wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

# MacOS
wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

# Conda installer (Linux and macOS)
conda install -c conda-forge -c bioconda foldseek
```
Other precompiled binaries for ARM64 amd SSE2 are available at [https://mmseqs.com/foldseek](https://mmseqs.com/foldseek).

## Memory requirments 
For optimal software performance, consider three options based on your RAM and search requirements:

1. **With Cα info (default).** 
   Use this formula to calculate RAM - `(6 bytes Cα + 1 3Di byte + 1 AA byte) * (database residues)`. The 54M AFDB50 entries require 151GB.

2. **Without Cα info.** 
   By disabling `--sort-by-structure-bits 0`, RAM requirement reduces to 35GB. However, this alters hit rankings and final scores but not E-values. Structure bits are mostly relevant for hit ranking for E-value > 10^-1.

3. **Single query searches.** 
   Use the `--prefilter-mode 1`, which isn't memory-limited and computes all ungapped alignments. This option optimally utilizes foldseek's multithreading capabilities for single queries.

## Tutorial Video
We presented a Foldseek tutorial at the SBGrid where we demonstrate the webserver and command line interface of foldseek. 
Check it out [here](https://www.youtube.com/watch?v=k5Rbi22TtOA).

<a href="https://www.youtube.com/watch?v=k5Rbi22TtOA"><img src="https://img.shields.io/youtube/views/k5Rbi22TtOA?style=social"></a>.

## Documentation
Many of Foldseek's modules (subprograms) rely on MMseqs2. For more information about these modules, refer to the [MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki). For documentation specific to Foldseek, checkout the Foldseek wiki [here](https://github.com/steineggerlab/foldseek/wiki).

## Quick start

### Search
The `easy-search` module allows to search single or multiple query structures, formatted in PDB/mmCIF format (flat or gzipped), against a target database, folder or single protein structures. In default it outputs the alignment information as a [tab-separated file](#tab-separated) but we support also [Superposed Cα PDBs](#superpositioned-cα-only-pdb-files) or a [HTML](#interactive-html) output.

    foldseek easy-search example/d1asha_ example/ aln tmpFolder
    
#### Output Search
##### Tab-separated
  
The default fields are containing the following fields: `query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits` but they can be customized with the `--format-output` option e.g. `--format-output "query,target,qaln,taln"` returns the query and target accession and the pairwise alignments in tab separated format. You can choose many different output columns.

| Code | Description |
| --- | --- |
|query | Query sequence identifier |
|target | Target sequence identifier |
|qca        | Calpha coordinates of the query |
|tca        | Calpha coordinates of the target |
|alntmscore | TM-score of the alignment | 
|qtmscore   | TM-score normalized by the query length |
|ttmscore   | TM-score normalized by the target length |
|u          | Rotation matrix (computed to by TM-score) |
|t          | Translation vector (computed to by TM-score) |
|lddt       | Average LDDT of the alignment |
|lddtfull   | LDDT per aligned position |
|prob       | Estimated probability for query and target to be homologous (e.g. being within the same SCOPe superfamily) |

Check out the [MMseqs2 documentation for more format output codes](https://github.com/soedinglab/MMseqs2/wiki#custom-alignment-format-with-convertalis).

##### Superpositioned Cα only PDB files
Foldseek's `--format-mode 5` generates PDB files with all Cα atoms superimposed based on the aligned coordinates on to the query structure. 
For each pairwise alignment it will write a single PDB files, so be carefull when using this options for large searches. 

##### Interactive HTML
Foldseek can locally generate a search result HTML similiar to the [webserver](https://search.foldseek.com) by specifying the format mode `--format-mode 3`

```
foldseek easy-search example/d1asha_ example/ result.html tmp --format-mode 3
```

<p align="center"><img src="./.github/results.png" height="400"/></p>

#### Important search parameters

| Option            | Category        | Description                                                                                               |
|-------------------|-----------------|-----------------------------------------------------------------------------------------------------------|
| -s              | Sensitivity     | Adjust sensitivity to speed trade-off; lower is faster, higher more sensitive (fast: 7.5, default: 9.5)   |
| --exhaustive-search | Sensitivity | Skips prefilter and performs an all-vs-all alignment (more sensitive but much slower)                     |
| --max-seqs      | Sensitivity     | Adjust the amount of prefilter handed to alignment; increasing it can lead to more hits (default: 1000)   |
| -e              | Sensitivity     | List matches below this E-value (range 0.0-inf, default: 0.001); increasing it reports more distant structures |
| --alignment-type| Alignment       | 0: 3Di Gotoh-Smith-Waterman (local, not recommended), 1: TMalign (global, slow), 2: 3Di+AA Gotoh-Smith-Waterman (local, default) |
| -c              | Alignment  | List matches above this fraction of aligned (covered) residues (see --cov-mode) (default: 0.0); higher coverage = more global alignment |
| --cov-mode      | Alignment  | 0: coverage of query and target, 1: coverage of target, 2: coverage of query                               |

#### Alignment Mode
In default Foldseek uses its local 3Di+AA strutural alignment but it also supports to realign hits using the global TMalign as well as rescoring alignments using TMscore. 

    foldseek easy-search example/d1asha_ example/ aln tmp --alignment-type 1

In case of the alignment type (`--alignment-type 1`) tmalign, we sort the results by the TMscore normalized by query length. We write the TMscore into the e-value=(qTMscore+tTMscore)/2 as well as into the score(=qTMscore*100) field. All output fields (like pident, fident, and alnlen) are calculated from the TMalign alignment.

### Databases 
The `databases` command downloads pre-generated databases like PDB or AlphaFoldDB.
    
    # pdb  
    foldseek databases PDB100 pdb tmp 
    # alphafold db
    foldseek databases Alphafold/Proteome afdb tmp 

We currently support the following databases: 
```
  Name                   	Type     	Taxonomy	Url
- Alphafold/UniProt   	Aminoacid	     yes	https://alphafold.ebi.ac.uk/
- Alphafold/UniProt50 	Aminoacid	     yes	https://alphafold.ebi.ac.uk/
- Alphafold/Proteome  	Aminoacid	     yes	https://alphafold.ebi.ac.uk/
- Alphafold/Swiss-Prot	Aminoacid	     yes	https://alphafold.ebi.ac.uk/
- ESMAtlas30          	Aminoacid	       -	https://esmatlas.com
- PDB                 	Aminoacid	     yes	https://www.rcsb.org
```

#### Create custom databases and indexes
The target database can be pre-processed by `createdb`. This make sense if searched multiple times. 
 
    foldseek createdb example/ targetDB
    foldseek createindex targetDB tmp  #OPTIONAL generates and stores the index on disk
    foldseek easy-search example/d1asha_ targetDB aln.m8 tmpFolder

### Cluster
The `easy-cluster` algorithm is designed for structural clustering by assigning structures to a representative protein using structural alignment. It accepts input in either PDB or mmCIF format, with support for both flat and gzipped files. By default, easy-cluster generates three output files with the following prefixes: (1) `_clu.tsv`, (2) `_repseq.fasta`, and (3) `_allseq.fasta`. The first file (1) is a [tab-separated](#tab-separated-cluster) file describing the mapping from representative to member, while the second file (2) contains only [representative sequences](#representative-fasta), and the third file (3) includes all [cluster member sequences](#all-member-fasta).

    foldseek easy-cluster example/ res tmp -c 0.9 
    
#### Output Cluster
##### Tab-separated cluster
The provided format represents protein structure clustering in a tab-separated, two-column layout (representative and member). Each line denotes a cluster-representative and cluster-member relationship, signifying that the member shares significant structural similarity with the representative, and thus belongs to the same cluster.
```
Q0KJ32	Q0KJ32
Q0KJ32	C0W539
Q0KJ32	D6KVP9
E3HQM9	E3HQM9
E3HQM9	F0YHT8
```

##### Representative fasta
The `_repseq.fasta` contains all representative protein sequences of the clustering.
```
>Q0KJ32
MAGA....R
>E3HQM9
MCAT...Q
```

##### All member fasta
In `_allseq.fasta` file all sequences of the cluster are present. A new cluster is marked by two identical name lines of the representative sequence, where the first line stands for the cluster and the second is the name line of the first cluster sequence. It is followed by the fasta formatted sequences of all its members.

```
>Q0KJ32	
>Q0KJ32
MAGA....R
>C0W539
MVGA....R
>D6KVP9
MVGA....R
>D1Y890
MVGV....R
>E3HQM9	
>E3HQM9
MCAT...Q
>Q223C0
MCAR...Q
```

#### Important cluster parameters

| Option            | Category        | Description                                                                                               |
|-------------------|-----------------|-----------------------------------------------------------------------------------------------------------|
| -e              | Sensitivity     | List matches below this E-value (range 0.0-inf, default: 0.001); increasing it reports more distant structures |
| --alignment-type| Alignment       | 0: 3Di Gotoh-Smith-Waterman (local, not recommended), 1: TMalign (global, slow), 2: 3Di+AA Gotoh-Smith-Waterman (local, default) |
| -c              | Alignment  | List matches above this fraction of aligned (covered) residues (see --cov-mode) (default: 0.0); higher coverage = more global alignment |
| --cov-mode      | Alignment  | 0: coverage of query and target, 1: coverage of target, 2: coverage of query                               |
| --min-seq-id      | Alignment  | the minimum sequence identity to be clustered                               |


## Main Modules
- `easy-search`       fast protein structure search  
- `easy-cluster`      fast protein structure clustering  
- `createdb`          create a database from protein structures (PDB,mmCIF, mmJSON)
- `databases`         download pre-assembled databases

## Examples
### Rescore aligments using TMscore
Easiest way to get the alignment TMscore normalized by min(alnLen,qLen,targetLen) as well as a rotation matrix is through the following command:
```
foldseek easy-search example/ example/ aln tmp --format-output query,target,alntmscore,u,t
```

Alternative, it is possible to compute TMscores for the kind of alignment output (e.g. 3Di/AA) using the following commands: 
```
foldseek createdb example/ targetDB
foldseek createdb example/ queryDB
foldseek search queryDB targetDB aln tmpFolder -a
foldseek aln2tmscore queryDB targetDB aln aln_tmscore
foldseek createtsv queryDB targetDB aln_tmscore aln_tmscore.tsv
```

Output format `aln_tmscore.tsv`: query and target identifier, TMscore, translation(3) and rotation vector=(3x3)

### Cluster search results 
The following command aligns the input structures all-against-all and keeps only alignments with 80% of the sequence covered by the alignment (-c 0.8) (read more about alignment coverage [here](https://github.com/soedinglab/MMseqs2/wiki#how-to-set-the-right-alignment-coverage-to-cluster)). It then clusters the results using greedy set cover algorithm. The clustering mode can be adjusted using --cluster-mode, read more [here](https://github.com/soedinglab/MMseqs2/wiki#clustering-modes). The clustering output format is described [here](https://github.com/soedinglab/MMseqs2/wiki#cluster-tsv-format).

```
foldseek createdb example/ db
foldseek search db db aln tmpFolder -c 0.8 
foldseek clust db aln clu
foldseek createtsv db db clu clu.tsv
```

### Query centered multiple sequence alignment 
Foldseek can generate a3m based multiple sequence alignments using the following commands. 
a3m can be converted to fasta format using [reformat.pl](https://raw.githubusercontent.com/soedinglab/hh-suite/master/scripts/reformat.pl) (`reformat.pl in.a3m out.fas`).
```
foldseek createdb example/ targetDB
foldseek createdb example/ queryDB
foldseek search queryDB targetDB aln tmpFolder -a
foldseek result2msa queryDB targetDB aln msa --msa-format-mode 6
foldseek unpackdb msa msa_output --unpack-suffix a3m --unpack-name-mode 0
```
