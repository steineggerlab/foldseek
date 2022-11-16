# Foldcomp
<p align="center">
<img src="https://raw.githubusercontent.com/steineggerlab/foldcomp/master/.github/img/foldcomp_strong_marv.png" max-height="300px" height="300" display="block" margin-left="auto" margin-right="auto" display="block"/>
</p>
Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes per residue, thus an averaged-sized protein of 350 residues requires ~6kb. 

Foldcomp efficient compressed format stores protein structures requiring only 13 bytes per residue, which reduces the required storage space by an order of magnitude compared to saving 3D coordinates directly. We achieve this reduction by encoding the torsion angles of the backbone as well as the side-chain angles in a compact binary file format (FCZ).

> Foldcomp currently only supports compression of single chain PDB files
<br clear="right"/>

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/steineggerlab/foldcomp/master/.github/img/format_benchmark_dark.png">
  <img src="https://raw.githubusercontent.com/steineggerlab/foldcomp/master/.github/img/format_benchmark_light.png" alt="Left panel: Foldcomp data format, saving amino acid residue in 13 byte. Top right panel:  Foldcomp decompression is as fast as gzip. Bottom right panel: Foldcomp compression ratio is higher than pulchra and gzip." max-width="720px" max-height="400px" width="auto" height="auto">
</picture>
</p>

## Usage

### Installing Foldcomp

```
# Install Foldcomp Python package
pip install foldcomp

# Download static binaries for Linux
wget https://mmseqs.com/foldcomp/foldcomp-linux-x86_64.tar.gz

# Download static binaries for Linux (ARM64)
wget https://mmseqs.com/foldcomp/foldcomp-linux-arm64.tar.gz

# Download binary for macOS
wget https://mmseqs.com/foldcomp/foldcomp-macos-universal.tar.gz
```

### Executable
```
# Compression
foldcomp compress <pdb_file|cif_file> [<fcz_file>]
foldcomp compress [-t number] <pdb_dir|cif_dir> [<fcz_dir>]

# Decompression
foldcomp decompress <fcz_file> [<pdb_file>]
foldcomp decompress [-t number] <fcz_dir> [<pdb_dir>]

# Extraction of sequence or pLDDT
foldcomp extract [--plddt|--fasta] <fcz_file> [<txt_file|fasta_file>]
foldcomp extract [--plddt|--fasta] [-t number] <fcz_dir|tar> [<output_dir>]

# Check
foldcomp check <fcz_file>
foldcomp check [-t number] <fcz_dir|tar>

# RMSD
foldcomp rmsd <pdb1|cif1> <pdb2|cif2>

# Options
 -h, --help           print this help message
 -t, --threads        threads for (de)compression of folders/tar files [default=1]
 -a, --alt            use alternative atom order [default=false]
 -b, --break          interval size to save absolute atom coordinates [default=25]
 -z, --tar            save as tar file [default=false]
 --plddt              extract pLDDT score (only for extraction mode)
 --fasta              extract amino acid sequence (only for extraction mode)
 --no-merge           do not merge output files (only for extraction mode)
```

### Python API

You can find more in-depth examples of using Foldcomp's Python interface in the example notebook:
<a href="https://colab.research.google.com/github/steineggerlab/foldcomp/blob/master/foldcomp-py-examples.ipynb" target="_blank" rel="noopener"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

```py
import foldcomp
# 01. Handling a FCZ file
# Open a fcz file
with open("test/compressed.fcz", "rb") as fcz:
  fcz_binary = fcz.read()

  # Decompress
  (name, pdb) = foldcomp.decompress(fcz_binary) # pdb_out[0]: file name, pdb_out[1]: pdb binary string

  # Save to a pdb file
  with open(name, "w") as pdb_file:
    pdb_file.write(pdb)

# 02. Iterate over a database of FCZ files
# Open a foldcomp database
ids = ["d1asha_", "d1it2a_"]
with foldcomp.open("test/example_db", ids=ids) as db:
  # Iterate through database
  for (name, pdb) in db:
      # save entries as seperate pdb files
      with open(name + ".pdb", "w") as pdb_file:
        pdb_file.write(pdb)
```

## Contributor
<a href="https://github.com/steineggerlab/foldcomp/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/foldcomp" />
</a>

