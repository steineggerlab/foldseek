# PDB server
SERVER="ftp.pdbj.org" # PDBj
CIF_SUBDIR="::ftp_data/structures/divided/mmCIF/" # subdirectory

# foldseek executable
foldseek="./foldseek/bin/foldseek"

# directories
SCRATCH="/scratch_dir_that_will_be_deleted_rmrf_after_completion"
MIRRORDIR="/parent_directory_to_download_all_pdb_files"

# pdb download log file
LOGFILE="file_to_record_logs_for_downloading_all_pdb_files"

# cloud storage upload endpoint
cloudflare_endpoint_url="endpoint_url_for_r2_cloudflarestorage_with_bucket_name=foldseek"
