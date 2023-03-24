# PDB server
SERVER="ftp.pdbj.org" # PDBj
CIF_SUBDIR="::ftp_data/structures/divided/mmCIF/" # subdirectory

# foldseek executable
foldseek="./foldseek/bin/foldseek"

# directories
MIRRORDIR="/parent_directory_to_download_all_pdb_files"
PAST_DB_FOLDER="/directory_to_store_previously_generated_foldseek_databases"
SYMLINK_FOLDER="/directory_for_symlinks_to_pdb_files"
# log file
LOGFILE="file_to_record_logs_for_downloading_all_pdb_files"

cloudflare_endpoint_url="endpoint_url_for_r2_cloudflarestorage_with_bucket_name=foldseek"
