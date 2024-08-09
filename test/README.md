# foldseek ProstT5 integration
For compiling foldseek binaries use the following conda environment
```
conda create -n fs-prostt5 -c conda-forge compilers cmake cuda rust openblas
```

To compile the foldseek binries run this command on a gpu server (e.g. devbox001)
```bash
srun -c 32 -t 20-0 --pty /bin/bash
# compiling with CUDA did not work, thus disable CUDA for now
# cmake -DCMAKE_INSTALL_PREFIX=. -DWITH_CUDA=1 -DENABLE_CUDA=1 -DCUDAToolkit_ROOT=/home/mpjw/miniforge3/envs/fs-prostt5 ..

# run this in foldseek repo root
mkdir -p build-test
cd build-test
cmake -DCMAKE_INSTALL_PREFIX=. -DENABLE_CUDA=0 -DIGNORE_RUST_VERSION=1 -DCMAKE_BUILD_TYPE=debug ..

make -j 32

ln -s ../test/ .
cd src/
ln -s ../../src/strucclustutils/structcreatedb.cpp
cd ..
```

## ProstT5 long sequence splitting 
Testing the splitting of long sequences functionality is designed for a max split length of 500.
Use the long sequences in `data/long.500.aa.fasta` and the pre split sequences in `data/short.500.aa.fasta`.
For each of this data create a foldseek database, extract the 3Di fasta data and compare it to each other (e.g. using [this python script](https://github.com/mpjw/ProstT5/blob/auto_split_long_seq/test/print_3Di_diff.py)).
If splitting works correctly 3Di sequences should observe 100% sequence identity.
```bash
# create folders for test databases
mkdir -p test/dbs/

# create foldseek db based on ProstT5 predicted 3Di
cd build-test
src/foldseek createdb test/data/long.500.aa.fasta test/dbs/long_500 --prostt5-model /home/sukhwan/foldseek_ctranslate/foldseek/weights/model --threads 32
src/foldseek lndb test/dbs/long_500_h test/dbs/long_500_ss_h
src/foldseek convert2fasta test/dbs/long_500_ss test/dbs/long.500.3Di.fasta

# create foldseek db with no splitting (all sequences shorter 500)
src/foldseek createdb test/data/short.500.aa.fasta test/dbs/short_500 --prostt5-model /home/sukhwan/foldseek_ctranslate/foldseek/weights/model --threads 32
src/foldseek lndb test/dbs/short_500_h test/dbs/short_500_ss_h
src/foldseek convert2fasta test/dbs/short_500_ss test/dbs/short.500.3Di.fasta
```
