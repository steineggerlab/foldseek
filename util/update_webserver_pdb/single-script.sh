#!/bin/bash
# Exit when any command fails
set -e
set -x

# Keep track of last executed command & echo error message before exit
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

source ./config.sh

# Check if foldseek is installed
if command -v $foldseek >/dev/null 2>&1; then
	echo "foldseek found"
else
	echo "foldseek not found"
	exit 1
fi

# Check if aws cli is installed
if command -v aws >/dev/null 2>&1; then
	echo "aws found"
else
	echo "aws not found"
	exit 1
fi

## 1) Restore backup of last database
lastDB=$(find pdb100.tar.gz)
thisDate=$(cat pdb_download_date)
mv $lastDB $lastDB$thisDate 
mv $lastDB$thisDate $PAST_DB_FOLDER
rm pdb*

## 2) Download mmCIF files from ftp.pdbj.org to PDB/FTP/
# rsync PDB files with accession number in string uniqueProteins, each separated by NEWLINE
date -u +%Y-%m-%d > pdb_download_date
rsync -av -zP ${SERVER}${CIF_SUBDIR} $MIRRORDIR > $LOGFILE 2>/dev/null

# store absolute paths of all files nested in PDB/FTP directory
find $MIRRORDIR -type f > full_filepaths.txt

## 3) Download ID Mapping file from ftp.uniprot.org to Taxonomy
rm -f zipped_idmap.gz idmap.txt
aria2c -x 8 --log=idmap_log --out=zipped_idmap.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
gzip -c -d zipped_idmap.gz > idmap.txt

## 4)  Process idmap.txt into two-column format (PDBID "\t" TaxonID)
awk -F"\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="BAD" } 1' OFS="\t" idmap.txt | awk -F '\t' '{
if ($6 != "BAD" && $13 != "BAD")
	print $6 "\t" $13
else if ($6 != "BAD" && $13 = "BAD")
	print $6 "\t" "0"
}' > idmap_columns_6_and_13

awk -F"\t" 'split($1,f,/; /)>1{for (i=1; i in f; i++) {$1=f[i]; print} next } 1' idmap_columns_6_and_13 > temp.txt
cat temp.txt | awk '{print tolower(substr($0,1,4)) substr($0,5,length($0)-4)}' | awk -F' ' 'BEGIN{OFS=FS} $2 == "" {$2 = "0"} 1' | sed 's/:/.cif.gz_/' > id_mapping_file 
echo "$(sort -u id_mapping_file)" > id_mapping_file

rm idmap_columns_6_and_13 temp.txt

# Create symlinks to each pdb file
while read p; do
	fileName=$(echo $p | awk -F"/" '{print $NF}')
	if ! test -f "${SYMLINK_FOLDER}/$fileName"; then
		ln -s $p ./symlinks
	fi
done < full_filepaths.txt

## 5) run foldseek createdb, cluster workflows
# Createdb
$foldseek createdb ./symlinks targetDB --chain-name-mode 1

# Get rid of .cif.gz from foldseek search results
$foldseek prefixid targetDB_h targetDB_h.tsv --tsv
awk '{gsub(".cif.gz",""); print $0}' targetDB_h.tsv > targetDB_h_new.tsv
$foldseek tsv2db targetDB_h_new.tsv targetDB_h  --output-dbtype 12

# Clustering: Create 'targetDB_clu100' by removing redundancies from 'targetDB'
$foldseek cluster targetDB targetDB_clu100 tmp -c 1.0 --min-seq-id 1.0 --cov-mode 0
# Create sub-database 'pdb', so we don't have to remove redundancies everytime we execute 'foldseek search'
$foldseek createsubdb targetDB_clu100 targetDB pdb
$foldseek createsubdb targetDB_clu100 targetDB_h pdb_h
$foldseek createsubdb targetDB_clu100 targetDB_ca pdb_ca
$foldseek createsubdb targetDB_clu100 targetDB_ss pdb_ss

## 6) Assign Taxonomy
$foldseek createtaxdb pdb tmp --tax-mapping-file id_mapping_file

# Append taxonomy ID = 0 to protein entries which are present in lookup file but missing in mapping file
TEMP_MAPPING="appendToMapping.txt"

MAX_ID=$(wc -l pdb.lookup | awk '{ print $1 }')
MAX_ID=$MAX_ID+1
LAST_OBSERVED_ID=-1

rm -f $TEMP_MAPPING
touch $TEMP_MAPPING

while read -r line
do
	THIS_ID=$(cut -f1 <<< "$line")
	for (( i=$LAST_OBSERVED_ID+1; i<$THIS_ID; i++ ))
	do
		echo "$i	0" >> $TEMP_MAPPING
	done
	LAST_OBSERVED_ID=$THIS_ID
done < pdb_mapping

cat $TEMP_MAPPING >> pdb_mapping
sort -k1 -n pdb_mapping -o pdb_mapping
rm $TEMP_MAPPING

# Convert softlinks to real files
SOFTLINKS=("pdb.lookup" "pdb.source" "pdb_ss_h" "pdb_ss_h.dbtype" "pdb_ss_h.index")
backup="-backup"

for softlink in ${SOFTLINKS[@]}
do
	if [ -e $softlink ]
	then
		cp $softlink $softlink$backup
		rm $softlink
		mv $softlink$backup $softlink
	fi
done

# Remove redundant indices from lookup, mapping files
DB="pdb"
indices="indices.txt"
lookupExtension=".lookup"
mappingExtension="_mapping"

awk -F"\t" '{print $1}' $DB.index > $indices

INFLATED=("$DB$lookupExtension" "$DB$mappingExtension")
extract="-extract"

for file in ${INFLATED[@]}
do
	echo $file
	awk -F"\t" 'NR==FNR{a[$1]; next} $1 in a' $indices $file > $file$extract
	rm $file
	mv $file$extract $file
done
rm $indices

# Remove useless database
rm targetDB*

## 7) Three example runs - sanity check
# Check if top result for 7p03 is 7p03, 3ojp is 3ojp, 9rub is 9rub
# If a redundant entry is added to PDB (ex: xyzw which is identical to 7p03),
# topID may become xyzw => sanity check failing although there's actually nothing wrong.
# If such an event occurs in the future, we shall simply replace p0/7p03 in example_runs with yz/xyzw 
declare -a example_runs=("p0/7p03" "oj/3ojp" "ru/9rub")
for id in "${example_runs[@]}"
do
	$foldseek easy-search ${MIRRORDIR}/${id}.cif.gz pdb aln.m8 tmp
	topResult=$(head -n 1 aln.m8 | awk -F"\t" '{print $2}')
	topID=${topResult:0:4} 
	readID=${id:3:4}
	if [ "$topID" != "$readID" ]; then
		echo "Failed sanity check with ${id}"
		exit 1
	fi
done

## 8) Upload to R2 storage
# Formatting: prepare upload to cloudflare
dbname="pdb100"
timestamp=$(cat pdb_download_date)
foldseekVersion=$($foldseek | grep "foldseek Version: " | sed 's/foldseek Version: //')
zip_filename=${dbname}.tar.gz
md5_file="pdb100.version"

md5sum pdb pdb_ca pdb_ca.dbtype pdb_ca.index pdb_h pdb_h.dbtype pdb_h.index pdb_mapping pdb_ss pdb_ss.dbtype pdb_ss.index pdb_taxonomy pdb.dbtype pdb.index pdb.lookup > pdb.md5sum
tar czf $zip_filename pdb pdb_ca pdb_ca.dbtype pdb_ca.index pdb_h pdb_h.dbtype pdb_h.index pdb_mapping pdb_ss pdb_ss.dbtype pdb_ss.index pdb_taxonomy pdb.dbtype pdb.index pdb.lookup pdb.md5sum
md5sum $zip_filename > $md5_file 
stamp="${timestamp}\tPDB_DATE\n${foldseekVersion}\tFOLDSEEK_COMMIT"
echo -e $stamp >> $md5_file 

# Make sure to not hit R2 upload rate limits
aws configure set default.s3.max_concurrent_requests 2
aws configure set default.s3.multipart_threshold 128MB
aws configure set default.s3.multipart_chunksize 128MB

# Upload zipped file and md5 file to cloudflare r2 storage
aws s3 cp $zip_filename s3://foldseek/$zip_filename --endpoint-url $cloudflare_endpoint_url
aws s3 cp $md5_file s3://foldseek/$md5_file --endpoint-url $cloudflare_endpoint_url
