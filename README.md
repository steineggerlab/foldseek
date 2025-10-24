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

## Publications

[Hyunbin Kim, Milot Mirdita, Martin Steinegger, Foldcomp: a library and format for compressing and indexing large protein structure sets, Bioinformatics, 2023;, btad153,](https://doi.org/10.1093/bioinformatics/btad153)

## Presentation Video

We presented Foldcomp at ISMB/ECCB2023. Check it out:

<a href="https://www.youtube.com/watch?v=aFtqH0VqE7w" target="_blank">
  <img src="https://raw.githubusercontent.com/steineggerlab/foldcomp/master/.github/img/ismb_thumbnail.png" alt="Foldcomp presented at ISMB/ECCB2023" max-width="720px" max-height="400px" width="auto" height="auto">
</a>

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

# Download binary for Windows (x64)
wget https://mmseqs.com/foldcomp/foldcomp-windows-x64.zip
```

### Executable
```
# Compression
foldcomp compress <pdb|cif> [<fcz>]
foldcomp compress [-t number] <dir|tar(.gz)> [<dir|tar|db>]

# Decompression
foldcomp decompress <fcz|tar> [<pdb>]
foldcomp decompress [-t number] <dir|tar(.gz)|db> [<dir|tar>]

# Decompressing a subset of Foldcomp database
foldcomp decompress [-t number] --id-list <idlist.txt> <db> [<dir|tar>]

# Extraction of sequence or pLDDT
foldcomp extract [--plddt|--amino-acid] <fcz> [<fasta>]
foldcomp extract [--plddt|--amino-acid] [-t number] <dir|tar(.gz)|db> [<fasta_out>]

# Check
foldcomp check <fcz>
foldcomp check [-t number] <dir|tar(.gz)|db>

# RMSD
foldcomp rmsd <pdb|cif> <pdb|cif>

# Options
 -h, --help               print this help message
 -v, --version            print version
 -t, --threads            threads for (de)compression of folders/tar files [default=1]
 -r, --recursive          recursively look for files in directory [default=0]
 -f, --file               input is a list of files [default=0]
 -a, --alt                use alternative atom order [default=false]
 -b, --break              interval size to save absolute atom coordinates [default=25]
 -z, --tar                save as tar file [default=false]
 -d, --db                 save as database [default=false]
 -y, --overwrite          overwrite existing files [default=false]
 -l, --id-list            a file of id list to be processed (only for database input)
 --skip-discontinuous     skip PDB with with discontinuous residues (only batch compression)
 --check                  check FCZ before and skip entries with error (only for batch decompression)
 --plddt                  extract pLDDT score (only for extraction mode)
 -p, --plddt-digits       extract pLDDT score with specified number of digits (only for extraction mode)
                          - 1: single digit (fasta-like format), 2: 2-digit(00-99; tsv), 3: 3-digit, 4: 4-digit (max)
 --fasta, --amino-acid    extract amino acid sequence (only for extraction mode)
 --no-merge               do not merge output files (only for extraction mode)
 --use-title              use TITLE as the output file name (only for extraction mode)
 --time                   measure time for compression/decompression
```

### Downloading Databases
We offer prebuilt databases for multiple large sets of predicted protein structures and a Python helper to download the database files.

You can download the AlphaFoldDB Swiss-Prot with the following command:
```
python -c "import foldcomp; foldcomp.setup('afdb_swissprot_v4');
```

Currently we offer the following databases:
* [ESMAtlas](https://esmatlas.com/) full (v0 + v2023_02): `foldcomp.setup('esmatlas')`
* ESMAtlas v2023_02: `foldcomp.setup('esmatlas_v2023_02')`
* ESMAtlas high-quality: `foldcomp.setup('highquality_clust30')`

  **Note:** We skipped all structures with discontinous residues or other issues.

* [AlphaFoldDB Uniprot](https://alphafold.ebi.ac.uk/): `foldcomp.setup('afdb_uniprot_v4')`
* AlphaFoldDB Swiss-Prot: `foldcomp.setup('afdb_swissprot_v4')`
* AlphaFoldDB Model Organisms: `foldcomp.setup('h_sapiens')`
  * `a_thaliana`, `c_albicans`, `c_elegans`, `d_discoideum`, `d_melanogaster`, `d_rerio`, `e_coli`, `g_max`,
    `h_sapiens`, `m_jannaschii`, `m_musculus`, `o_sativa`, `r_norvegicus`, `s_cerevisiae`, `s_pombe`, `z_mays`
* [AlphaFoldDB Cluster Representatives](https://afdb-cluster.steineggerlab.workers.dev/): `foldcomp.setup('afdb_rep_v4')`
* AlphaFoldDB Cluster Representatives (Dark Clusters): `foldcomp.setup('afdb_rep_dark_v4')`

If you want other prebuilt datasets, please get in touch with us through our [GitHub issues](https://github.com/steineggerlab/foldcomp/issues).

If you have issues downloading the databases you can navigate directly to our [download server](https://foldcomp.steineggerlab.workers.dev/) and download the required files. E.g. `afdb_uniprot_v4`, `afdb_uniprot_v4.index`, `afdb_uniprot_v4.dbtype`, `afdb_uniprot_v4.lookup`, and optionally `afdb_uniprot_v4.source`.

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

  # Get data as dictionary
  data_dict = foldcomp.get_data(fcz_binary) # foldcomp.get_data(pdb) also works
  # Keys: phi, psi, omega, torsion_angles, residues, bond_angles, coordinates
  data_dict["phi"] # phi angles (C-N-CA-C)
  data_dict["psi"] # psi angles (N-CA-C-N)
  data_dict["omega"] # omega angles (CA-C-N-CA)
  data_dict["torsion_angles"] # torsion angles of the backbone as list (phi + psi + omega)
  data_dict["bond_angles"] # bond angles of the backbone as list
  data_dict["residues"] # amino acid residues as string
  data_dict["coordinates"] # coordinates of the backbone as list

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

## Subsetting Databases
If you are dealing with millions of entries, we recommend using `createsubdb` command
of [mmseqs2](https://mmseqs.com) to subset databases.
The following commands can be used to subset the AlphaFold Uniprot DB with given IDs.
```sh
# mmseqs createsubdb --subdb-mode 0 --id-mode 1 id_list.txt input_foldcomp_db output_foldcomp_db
mmseqs createsubdb --subdb-mode 0 --id-mode 1 id_list.txt afdb_uniprot_v4 afdb_subset
```
Please note that the IDs in afdb_uniprot_v4 are in the format `AF-A0A5S3Y9Q7-F1-model_v4` .

## Community Contributions
* [PyMOL Plugin for reading Foldcomp files](https://github.com/yakomaxa/load_fcz_PyMOL) by @yakomaxa

## Contributor
<a href="https://github.com/steineggerlab/foldcomp/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/foldcomp" />
</a>

