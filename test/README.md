# foldseek ProstT5 integration
This repository is part of the integration of [ProstT5](https://www.biorxiv.org/content/10.1101/2023.07.23.550085v1) into [foldseek](https://github.com/steineggerlab/foldseek).
More specifically, the split-wise prediction functionality is implemented here.
That is, splitting amino acid sequences longer than ProstT5 attention into smaller subsequences (maximum length of 6000 per default) and concatenating the 3Di representation sequences afterwards.

<p align="center"><img src="https://github.com/steineggerlab/foldseek/blob/master/.github/foldseek.png" height="250"/></p>

## Relevant Papers
[van Kempen M, Kim S, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, and Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)](https://www.nature.com/articles/s41587-023-01773-0)

[Heinzinger M, Weissenow K, Joaquin G S, Henkel A, Steinegger M and Rost B. ProstT5: Bilingual Language Model for Protein Sequence and Structure. bioRxiv, doi:10.1101/2023.07.23.550085 (2023)](https://www.biorxiv.org/content/10.1101/2023.07.23.550085v1)

## Compiling
For compiling foldseek binaries use the following conda environment
```bash
conda create -n cuda-dev -c conda-forge cmake cuda-nvcc rust openblas libcublas-dev libcublas-static "gcc<13" "gxx<13" libcurand-dev libcurand-static

# alternatively create from .yml file
conda env create -f cuda-dev.yml
```

To compile the foldseek binries run this command on a gpu server (e.g. devbox001)
```bash
# run this in foldseek repo root
mkdir -p build
cd build

cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. -DWITH_CUDA=1 -DENABLE_CUDA=1 -DCUDA_TOOLKIT_ROOT_DIR=$HOME/miniforge3/envs/cuda-dev/targets/x86_64-linux ..
# cmake -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_INSTALL_PREFIX=. -DENABLE_CUDA=0 -DIGNORE_RUST_VERSION=1 ..

srun -c 32 -t 20-0 --pty /bin/bash
make -j 32

# link test folder for debugging 
ln -s ../test/ .
cd src/
ln -s ../../src/strucclustutils/structcreatedb.cpp
cd ..
```

## Testing 
Data for testing the splitting of long sequences functionality of the ProstT5 integration in foldseek is provided in `test/data`.
Here, half of the sequences in `test/data/long.aa.fasta` are longer than ProstT5 attention.
Whereas `test/data/short.500.aa.fasta` and `test/data/short.6000.aa.fasta` contain the same sequences but pre-split with a maximum split length of `500` and `6000` respectively.
For testting, create a foldseek database for each split length pair, extract the 3Di fasta data and compare it (e.g. using [this python script](https://github.com/mpjw/ProstT5/blob/auto_split_long_seq/test/print_3Di_diff.py)).
If splitting works correctly 3Di sequences should observe 100% sequence identity.
```bash
# create folders for test databases
mkdir -p test/dbs/

# create foldseek db based on ProstT5 predicted 3Di
cd build-test
src/foldseek createdb test/data/long.500.aa.fasta test/dbs/long_500 --prostt5-model /home/sukhwan/foldseek_ctranslate/foldseek/weights/model --prostt5-split-length 500 --threads 32
src/foldseek lndb test/dbs/long_500_h test/dbs/long_500_ss_h
src/foldseek convert2fasta test/dbs/long_500_ss test/dbs/long.500.3Di.fasta

# create foldseek db from pre-splitted sequences
src/foldseek createdb test/data/short.500.aa.fasta test/dbs/short_500 --prostt5-model /home/sukhwan/foldseek_ctranslate/foldseek/weights/model --prostt5-split-length 0 --threads 32
src/foldseek lndb test/dbs/short_500_h test/dbs/short_500_ss_h
src/foldseek convert2fasta test/dbs/short_500_ss test/dbs/short.500.3Di.fasta

# create foldseek db based on ProstT5 predicted 3Di (ProstT5 implementation will do the splitting)
src/foldseek createdb test/data/long.6000.aa.fasta test/dbs/long_6000 --prostt5-model /home/sukhwan/foldseek_ctranslate/foldseek/weights/model --prostt5-split-length 6000 --threads 32
src/foldseek lndb test/dbs/long_6000_h test/dbs/long_6000_ss_h
src/foldseek convert2fasta test/dbs/long_6000_ss test/dbs/long.6000.3Di.fasta

# create foldseek db from pre-splitted sequences
src/foldseek createdb test/data/short.6000.aa.fasta test/dbs/short_6000 --prostt5-model /home/sukhwan/foldseek_ctranslate/foldseek/weights/model --prostt5-split-length 0 --threads 32
src/foldseek lndb test/dbs/short_6000_h test/dbs/short_6000_ss_h
src/foldseek convert2fasta test/dbs/short_6000_ss test/dbs/short.6000.3Di.fasta
```
