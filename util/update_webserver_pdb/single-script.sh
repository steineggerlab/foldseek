#!/bin/bash
# Exit when any command fails
set -e
set -x
set -u

last_command=""
current_command=""
# Keep track of last executed command & echo error message before exit
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "\"${last_command}\" command filed with exit code $?"' EXIT

source ./config.sh

# Check if foldseek is installed
if command -v $foldseek >/dev/null 2>&1; then
	echo "foldseek found"
else
	echo "foldseek not found"
	exit 1
fi
foldseekVersion=$($foldseek version)

# Check if aws cli is installed
if command -v aws >/dev/null 2>&1; then
	echo "aws found"
else
	echo "aws not found"
	exit 1
fi

# Download mmCIF files from ftp.pdbj.org to PDB/FTP/
# rsync PDB files with accession number in string uniqueProteins, each separated by NEWLINE
#rsync -av -zP ${SERVER}${CIF_SUBDIR} $MIRRORDIR > $LOGFILE 2>/dev/null
#date -u +%Y-%m-%d > pdb_download_date
timestamp=$(cat pdb_download_date)

mkdir -p ${SCRATCH}
pushd ${SCRATCH}
# Run foldseek createdb, cluster workflows
$foldseek createdb $MIRRORDIR pdb_seq --chain-name-mode 1 --write-mapping 1
# Get rid of .cif.gz from foldseek search results
$foldseek prefixid pdb_seq_h pdb_seq_h.tsv --tsv
sed -i 's|\.cif\.gz||g' pdb_seq_h.tsv 
sed -i 's|\.cif\.gz||g' pdb_seq.lookup 
sed -i 's|\.cif\.gz||g' pdb_seq.source
$foldseek tsv2db pdb_seq_h.tsv pdb_seq_h --output-dbtype 12
rm -f -- pdb_seq_h.tsv

# Clustering: Create 'targetDB_clu100' by removing redundancies from 'targetDB'
MMSEQS_FORCE_MERGE=1 $foldseek cluster pdb_seq pdb_clu tmp -c 0.95 --min-seq-id 1.0 --cov-mode 0
# Create sub-database 'pdb'
$foldseek createsubdb pdb_clu pdb_seq pdb --subdb-mode 1
rm pdb_h*
$foldseek lndb pdb_seq_h pdb_h
$foldseek createsubdb pdb_clu pdb_seq_ca pdb_ca --subdb-mode 1
$foldseek createsubdb pdb_clu pdb_seq_ss pdb_ss --subdb-mode 1
# Assign Taxonomy
$foldseek createtaxdb pdb_seq tmp
ln -sf -- pdb_seq.lookup pdb.lookup
ln -sf -- pdb_seq_mapping pdb_mapping
ln -sf -- pdb_seq_taxonomy pdb_taxonomy

# turn absolute symlinks to relative ones so we can make a correctly working tar
while read LINK TARGET; do
    rm -f -- $(basename $LINK)
    ln -sf -- $(basename $TARGET) $(basename $LINK)
done < <(find . -maxdepth 1 -type l -printf "%p\t%l\n")

# Three example runs - sanity check
# Check if top result for 7p03 is 7p03, 3ojp is 3ojp, 9rub is 9rub
# If a redundant entry is added to PDB (ex: xyzw which is identical to 7p03),
# topID may become xyzw => sanity check failing although there's actually nothing wrong.
# If such an event occurs in the future, we shall simply replace p0/7p03 in example_runs with yz/xyzw 
if false; then
declare -a example_runs=("p0/7p03" "oj/3ojp" "ru/9rub")
for id in "${example_runs[@]}"; do
	$foldseek easy-search ${MIRRORDIR}/${id}.cif.gz pdb aln.m8 tmp
	topResult=$(head -n 1 aln.m8 | awk -F"\t" '{print $2}')
	topID=${topResult:0:4} 
	readID=${id:3:4}
	if [ "$topID" != "$readID" ]; then
		echo "Failed sanity check with ${id}"
		exit 1
	fi
done
fi

## 8) Upload to R2 storage
# Formatting: prepare upload to cloudflare
md5sum pdb* > pdb.md5sum
tar --owner=0 --group=0 -czvf pdb100.tar.gz pdb* 
md5sum pdb100.tar.gz > pdb100.version
echo -e "${timestamp}\tPDB_DATE\n${foldseekVersion}\tFOLDSEEK_COMMIT" >> pdb100.version 
popd

cp -f -- ${SCRATCH}/pdb100.tar.gz .
cp -f -- ${SCRATCH}/pdb100.version .
rm -rf -- "${SCRATCH}"

# Make sure to not hit R2 upload rate limits
aws configure set default.s3.max_concurrent_requests 2
aws configure set default.s3.multipart_threshold 128MB
aws configure set default.s3.multipart_chunksize 128MB

# Upload zipped file and md5 file to cloudflare r2 storage
aws s3 cp pdb100.tar.gz s3://foldseek/pdb100.tar.gz --endpoint-url $cloudflare_endpoint_url
aws s3 cp pdb100.version s3://foldseek/pdb100.version --endpoint-url $cloudflare_endpoint_url
