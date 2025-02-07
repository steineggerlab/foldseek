#include "DownloadDatabase.h"
#include "Prefiltering.h"
#include "CommandDeclarations.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "structdatabases.sh.h"

LocalParameters& localPar = LocalParameters::getLocalInstance();
void updateValdiation() {
    DbValidator::allDb.push_back(LocalParameters::DBTYPE_CA_ALPHA);
    DbValidator::allDb.push_back(LocalParameters::DBTYPE_TMSCORE);
    DbValidator::allDbAndFlat.push_back(LocalParameters::DBTYPE_CA_ALPHA);
    DbValidator::allDbAndFlat.push_back(LocalParameters::DBTYPE_TMSCORE);
}
void (*validatorUpdate)(void) = updateValdiation;

std::vector<Command> foldseekCommands = {
        {"createdb",             structcreatedb,                &localPar.structurecreatedb,    COMMAND_MAIN,
                "Convert PDB/mmCIF/tar[.gz]/DB files or directory/TSV to a structure DB",
                "# Process multiple files\n"
                "foldseek createdb examples/1tim.pdb.gz examples/8tim.pdb.gz DB\n"
                "# Process a directory containing PDB|mmCIF[.gz]|tar[.gz]|DB recursively, only one directory can be given\n"
                "foldseek createdb examples/ DB\n"
                "# Process a TSV file with a list of PDB|mmCIF[.gz]|tar[.gz]|DB, only one TSV can be given\n"
                "foldseek createdb examples.tsv DB\n"
                "# Process a directory or tar file and filter based on file name\n"
                "# Note: --file-include and --file-exclude only apply to directory or tar input\n"
                "foldseek createdb examples/ DB --file-include \"pdb.gz$\"\n"
                "# Predict 3Di sequences from an amino acid FASTA file using ProstT5\n"
                "foldseek databases ProstT5 weights tmp\n"
                "foldseek createdb QUERY.fasta DB --prostt5-model weights\n\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:directory|.tsv>|<i:PDB|mmCIF[.gz]|tar[.gz]|DB> ... <i:PDB|mmCIF[.gz]|tar|DB> <o:sequenceDB>",
                CITATION_FOLDSEEK | CITATION_PROSTT5, {{"PDB|mmCIF[.gz]|stdin|tar[.gz]|DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC,
#ifdef HAVE_GCS
                                            &DbValidator::flatfileStdinGenericUri
#else
                                            &DbValidator::flatfileStdinAndGeneric
#endif
                                    },
                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"structureto3didescriptor",             structureto3didescriptor,            &localPar.structurecreatedb,    COMMAND_HIDDEN,
                "Convert PDB/mmCIF/tar[.gz] files to a db",
                "Convert PDB/mmCIF/tar[.gz] files to a db",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]> <o:3didescriptor>",
                CITATION_FOLDSEEK, {{"PDB|mmCIF[.gz]|stdin", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfileStdinAndGeneric },
                                           {"3didescriptor", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"createsubdb",          createstructsubdb,          &localPar.createstructsubdb,          COMMAND_SET,
                "Create a subset of a DB from list of DB keys",
                "# Create a new sequence, 3di, c-alpha DB with keys 1, 2 and 3\n"
                "foldseek createsubdb <(printf '1\n2\n3\n') sequenceDB oneTwoThreeDB\n\n"
                "# Create a new sequence, 3di, c-alpha database with representatives of clusterDB\n"
                "foldseek cluster sequenceDB clusterDB tmp\n"
                "foldseek createsubdb clusterDB sequenceDB representativesDB\n",
                "Milot Mirdita <milot@mirdita.de> & Sooyoung Cha <ellen2g77@gmail.com>",
                "<i:subsetFile|DB> <i:DB> <o:DB>",
                CITATION_FOLDSEEK|CITATION_MMSEQS2, {{"subsetFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDbAndFlat },
                                          {"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"easy-search",          easystructuresearch,           &localPar.easystructuresearchworkflow,   COMMAND_EASY,
                "Structual search",
                "# Search a single/multiple PDB file against a set of PDB files\n"
                "foldseek easy-search examples/d1asha_ examples/ result.m8 tmp\n"
                "# Format output differently\n"
                "foldseek easy-search examples/d1asha_ examples/ result.m8 tmp --format-output query,target,qstart,tstart,cigar\n"
                "# Align with TMalign (global)\n"
                "foldseek easy-search examples/d1asha_ examples/ result.m8 tmp --alignment-type 1\n"
                "# Skip prefilter and perform an exhaustive alignment (slower but more sensitive)\n"
                "foldseek easy-search examples/d1asha_ examples/ result.m8 tmp --exhaustive-search 1\n\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]>|<i:stdin> <i:targetFastaFile[.gz]>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_FOLDSEEK, {{"PDB|mmCIF[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &FoldSeekDbValidator::flatfileStdinAndFolder },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::flatfileAndFolder },
                                           {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-cluster",         easystructurecluster,          &localPar.easystructureclusterworkflow, COMMAND_EASY,
                "Slower, sensitive clustering",
                "foldseek easy-cluster examples/ result tmp\n"
                "# Cluster output\n"
                "#  - result_rep_seq.fasta: Representatives\n"
                "#  - result_all_seq.fasta: FASTA-like per cluster\n"
                "#  - result_cluster.tsv:   Adjacency list\n\n"
                "# Important parameter: --min-seq-id, --cov-mode and -c \n"
                "#                  --cov-mode \n"
                "#                  0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "#        -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]> <o:clusterPrefix> <tmpDir>",
                CITATION_FOLDSEEK, {{"PDB|mmCIF[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &FoldSeekDbValidator::flatfileStdinAndFolder },
                                           {"clusterPrefix", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"search",               structuresearch,               &localPar.structuresearchworkflow,       COMMAND_MAIN,
                "Sensitive homology search",
                "# Search multiple structures (mmCIF,PDB,tar) against structures.\n"
                "foldseek search queryDB targetDB resultDB tmp\n"
                "foldseek convertalis queryDB targetDB resultDB result.m8\n\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-rbh",                  structureeasyrbh,                  &localPar.easystructuresearchworkflow,       COMMAND_EASY,
                "Find reciprocal best hit",
                "# Assign reciprocal best hit\n"
                "mmseqs easy-rbh examples/QUERY.fasta examples/DB.fasta result tmp\n\n",
                "Eli Levy Karin & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryFastaFile1[.gz|.bz2]> <i:targetFastaFile[.gz|.bz2]>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_MMSEQS2,{{"fastaFile[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &FoldSeekDbValidator::flatfileStdinAndFolder },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::flatfileStdinAndFolder },
                                           {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"rbh",                  structurerbh,                  &localPar.structuresearchworkflow,       COMMAND_MAIN,
                "Reciprocal best hit search",
                NULL,
                "Eli Levy Karin & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"makepaddedseqdb",               makepaddeddb,              &localPar.makepaddeddb,              COMMAND_HIDDEN,
                "Generate a padded sequence, 3di, c-alpha DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr> & Sooyoung Cha <ellen2g77@gmail.com>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_GPU, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                          {"sequenceIndexDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"result2profile",       result2structprofile,       &localPar.result2structprofile,       COMMAND_PROFILE,
                "Compute profile DB from a result DB for both amino acid and 3di",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> & Sooyoung Cha <ellen2g77@gmail.com>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:profileDB>",
                CITATION_FOLDSEEK|CITATION_MMSEQS2,{{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"profileDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::profileDb }}},
        {"cluster",              structurecluster,   &localPar.structureclusterworkflow,      COMMAND_MAIN,
                "Slower, sensitive clustering",
                "# Cascaded clustering of FASTA file\n"
                "mmseqs cluster sequenceDB clusterDB tmp\n\n"
                "#                  --cov-mode \n"
                "# Sequence         0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "# Cutoff -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n"
                "# Cascaded clustering with reassignment\n"
                "# - Corrects criteria-violoations of cascaded merging\n"
                "# - Produces more clusters and is a bit slower\n"
                "mmseqs cluster sequenceDB clusterDB tmp --cluster-reassign\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr> & Lars von den Driesch",
                "<i:sequenceDB> <o:clusterDB> <tmpDir>",
                CITATION_FOLDSEEK|CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},



        {"tmalign",      tmalign,      &localPar.tmalign,      COMMAND_ALIGNMENT,
                "Compute tm-score ",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                           {"alnDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::alignmentDb }}},
        {"structurealign",      structurealign,      &localPar.structurealign,      COMMAND_ALIGNMENT,
                "Compute structural alignment using 3Di alphabet, amino acids and neighborhood information",
                NULL,
                "Charlotte Tumescheit <ch.tumescheit@gmail.com> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                           {"alnDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::alignmentDb }}},
        {"structurerescorediagonal",     structureungappedalign,       &localPar.structurerescorediagonal,      COMMAND_ALIGNMENT,
                "Compute sequence identity for diagonal",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"aln2tmscore", aln2tmscore,      &localPar.threadsandcompression,      COMMAND_ALIGNMENT,
                "Compute tmscore of an alignment database ",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:alnDB> <o:resultDB>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::sequenceDb },
                                           {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"tmDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::tmscore }}},
        {"samplemulambda", samplemulambda,      &localPar.samplemulambda,      COMMAND_EXPERT,
                "Sample mu and lambda from random shuffled sequences ",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:resultDB>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::sequenceDb },
                                           {"tmDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::genericDb }}},
        {"clust",                clust,                &localPar.clust,                COMMAND_CLUSTER,
                "Cluster result by Set-Cover/Connected-Component/Greedy-Incremental",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> & Lars von den Driesch & Maria Hauser",
                "<i:sequenceDB> <i:resultDB> <o:clusterDB>",
                CITATION_MMSEQS2|CITATION_MMSEQS1,{{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                           {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb }}},
        {"databases",            databases,            &localPar.databases,            COMMAND_DATABASE_CREATION,
                "List and download databases",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<name> <o:sequenceDB> <tmpDir>",
                CITATION_TAXONOMY|CITATION_FOLDSEEK, {{"selection", 0, DbType::ZERO_OR_ALL, &DbValidator::empty },
                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"tmpDir",     DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"createindex",          structureindex,       &localPar.createindex,          COMMAND_DATABASE_CREATION,
                "Store precomputed index on disk to reduce search overhead",
                "# Create protein sequence index\n"
                "mmseqs createindex sequenceDB tmp\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <tmpDir>",
                CITATION_SERVER | CITATION_FOLDSEEK,{{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"createclusearchdb", createclusearchdb,   &localPar.createclusearchdb,      COMMAND_DATABASE_CREATION,
                "Build a searchable cluster database allowing for faster searches",
                "# cluster database and build a searchable db\n"
                "foldseek cluster sequenceDB clusterDB tmp --min-seq-id 0.3\n"
                "foldseek createclusearchdb sequenceDB clusterDB clusterSearchDb\n"
                "foldseek search sequenceDB clusterSearchDb aln tmp --cluster-search 1\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <i:clusterDB> <o:sequenceDB>",
                CITATION_FOLDSEEK|CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"clusterDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                           {"clusterSearchDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"mmcreateindex",        createindex,          &localPar.createindex,          COMMAND_HIDDEN,
                NULL,
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <tmpDir>",
                CITATION_SERVER | CITATION_MMSEQS2,{{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"convertalis",          structureconvertalis,    &localPar.convertalignments,    COMMAND_FORMAT_CONVERSION,
                "Convert alignment DB to BLAST-tab, SAM or custom format",
                "# Create output in BLAST M8 format (12 columns):\n"
                "#  (1,2) identifiers for query and target sequences/profiles,\n"
                "#  (3) sequence identity, (4) alignment length, (5) number of mismatches,\n"
                "#  (6) number of gap openings, (7-8, 9-10) alignment start and end-position in query and in target,\n"
                "#  (11) E-value, and (12) bit score\n"
                "foldseek convertalis queryDB targetDB result.m8\n\n"
                "# Create a TSV containing pairwise alignments\n"
                "foldseek convertalis queryDB targetDB result.tsv --format-output query,target,qaln,taln\n\n"
                "# Annotate a alignment result with taxonomy information from targetDB\n"
                "foldseek convertalis queryDB targetDB result.tsv --format-output query,target,taxid,taxname,taxlineage\n\n"
                " Create SAM output\n"
                "foldseek convertalis queryDB targetDB result.sam --format-mode 1\n\n"
                "# Create a TSV containing which query file a result comes from\n"
                "foldseek createdb euk_queries.fasta bac_queries.fasta queryDB\n"
                "foldseek convertalis queryDB targetDB result.tsv --format-output qset,query,target\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDb> <i:targetDb> <i:alignmentDB> <o:alignmentFile>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"compressca",           compressca,             &localPar.compressca,            COMMAND_FORMAT_CONVERSION | COMMAND_EXPERT,
                "Create a new C-alpha DB with chosen compression encoding from a sequence DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:caDB>",
                CITATION_FOLDSEEK, {{"Db", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::sequenceDb },
                                           {"caDb", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::cadb }}},
        {"convert2pdb",          convert2pdb,             &localPar.convert2pdb,          COMMAND_FORMAT_CONVERSION,
                "Convert a foldseek structure db to a single multi model PDB file or a directory of PDB files",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:Db> <o:pdbFile|pdbDir>",
                CITATION_FOLDSEEK, {{"Db", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"pdbFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"scoremultimer", scoremultimer, &localPar.scoremultimer, COMMAND_ALIGNMENT,
                "Get multimer level alignments from alignmentDB",
                "# Get multimer level alignments (chain assignments and tm-scores) from alignmentDB.\n"
                "foldseek scoremultimer queryDB targetDB alignmentDB complexDB\n"
                "# simple tsv output format"
                "foldseek createmultimerreport queryDB targetDB complexDB result.tsv"
                "# output files with convertalis"
                "foldseek convertalis queryDB targetDB complexDB result.m8\n\n",
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:queryDb> <i:targetDb> <i:alignmentDB> <o:complexDB>",
                CITATION_FOLDSEEK_MULTIMER, {
                                           {"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::NEED_HEADER, &DbValidator::sequenceDb},
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::NEED_HEADER, &DbValidator::sequenceDb},
                                           {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb},
                                           {"complexDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb}
                                   }
        },
        {"scorecomplex", scoremultimer, &localPar.scoremultimer, COMMAND_HIDDEN,
                "", NULL, "", "", CITATION_FOLDSEEK_MULTIMER, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}
        },
        {"filtermultimer", filtermultimer, &localPar.filtermultimer, COMMAND_HIDDEN,
                "Filters multimers satisfying given coverage",
                "foldseek filtermultimer queryDB targetDB alignmentDB complexDB -c 0.8 --cov-mode 1\n",
                "Seongeun  Kim <seamustard52@gmail.com> & Sooyoung Cha <ellen2g77@gmail.com>",
                "<i:queryDB> <i:targetDB> <i:alignmentDB> <o:clustDB>",
                CITATION_FOLDSEEK_MULTIMER, {
                                           {"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"clustDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::clusterDb }
                                   }
        },
        {"multimercluster", multimercluster, &localPar.multimerclusterworkflow, COMMAND_MAIN, 
                "Multimer level cluster",
                "#Clustering of PDB DB\n"
                "foldseek multimercluster queryDB clusterDB tmp\n"
                "#                  --cov-mode \n"
                "# Sequence         0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "# Cutoff -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n",
                "Seongeun  Kim <seamustard52@gmail.com> & Sooyoung Cha <ellen2g77@gmail.com>",
                "<i:sequenceDB> <o:clusterDB> <tmpDir>",
                CITATION_FOLDSEEK_MULTIMER, {
                                        {"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                        {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::clusterDb },
                                        {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }
                }
        },     
        {"easy-multimercluster", easymultimercluster, &localPar.easymultimerclusterworkflow, COMMAND_EASY,
                "Multimer level cluster",
                "#Clustering of PDB files\n"
                "foldseek easy-multimercluster examples/ result tmp\n"
                "# Cluster output\n"
                "#  - result_rep_seq.fasta: Representatives\n"
                "#  - result_cluster.tsv:   Adjacency list\n\n"
                "# Important parameter: --cov-mode and -c \n"
                "#                  --cov-mode \n"
                "#                  0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "#        -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n",
                "Seongeun  Kim <seamustard52@gmail.com> & Sooyoung Cha <ellen2g77@gmail.com>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]> <o:clusterPrefix> <tmpDir>",
                CITATION_FOLDSEEK_MULTIMER, {
                                        {"PDB|mmCIF[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &FoldSeekDbValidator::flatfileStdinAndFolder},
                                        {"clusterPrefix", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                        {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }
                }
        },     
        {"multimersearch", multimersearch, &localPar.multimersearchworkflow, COMMAND_MAIN,
                "Multimer level search",
                "# Search a single/multiple PDB file against a set of PDB files and get multimer level alignments\n"
                "foldseek multimersearch queryDB targetDB result tmp\n"
                "# Format output differently\n"
                "foldseek multimersearch queryDB targetDB result tmp --format-output query,target,qstart,tstart,cigar\n"
                "# Align with TMalign (global)\n"
                "foldseek multimersearch queryDB targetDB result tmp --alignment-type 1\n"
                "# Skip prefilter and perform an exhaustive alignment (slower but more sensitive)\n"
                "foldseek multimersearch queryDB targetDB result tmp --exhaustive-search 1\n\n",
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_FOLDSEEK_MULTIMER, {
                                           {"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::NEED_HEADER, &DbValidator::sequenceDb},
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::NEED_HEADER, &DbValidator::sequenceDb},
                                           {"complexDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb},
                                           {"tempDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}
                                   }
        },
        {"complexsearch", multimersearch, &localPar.multimersearchworkflow, COMMAND_HIDDEN,
                "", NULL, "", "", CITATION_FOLDSEEK_MULTIMER, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}
        },
        {"easy-multimersearch", easymultimersearch, &localPar.easymultimersearchworkflow, COMMAND_EASY,
                "Multimer level search",
                "# Search a single/multiple PDB file against a set of PDB files and get multimer level alignments\n"
                "foldseek easy-multimersearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp\n"
                "# Format output differently\n"
                "foldseek easy-multimersearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp --format-output query,target,qstart,tstart,cigar\n"
                "# Align with TMalign (global)\n"
                "foldseek easy-multimersearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp --alignment-type 1\n"
                "# Skip prefilter and perform an exhaustive alignment (slower but more sensitive)\n"
                "foldseek easy-multimersearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp --exhaustive-search 1\n\n",
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]>|<i:stdin> <i:targetFastaFile[.gz]>|<i:targetDB> <o:outputFileName> <tmpDir>",
                CITATION_FOLDSEEK_MULTIMER, {
                                           {"PDB|mmCIF[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &FoldSeekDbValidator::flatfileStdinAndFolder},
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::flatfileAndFolder},
                                           {"outputFileName", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                           {"tempDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}
                                   }
        },
        {"easy-complexsearch", easymultimersearch, &localPar.easymultimersearchworkflow, COMMAND_HIDDEN,
                "", NULL, "", "", CITATION_FOLDSEEK_MULTIMER, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}
        },
        {"createmultimerreport", createmultimerreport, &localPar.createmultimerreport, COMMAND_FORMAT_CONVERSION,
                "Convert complexDB to tsv format",
                "# Create output in tsv format (9 columns):  qComplexName.c_str(), tComplexName.c_str(), qChainString.c_str(), tChainString.c_str(), qTMScore, tTMScore, u, t, assId\n"
                "#  (1,2) identifiers for query and target multimers,\n"
                "#  (3,4) chains of query multimer and target multimer,\n"
                "#  (5,6) tm score based on query and target residue length,\n"
                "#  (8,9) u and t,\n"
                "#  (9) assignment id\n"
                "foldseek createmultimerreport queryDB targetDB complexDB result.tsv\n",
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:queryDb> <i:targetDb> <i:complexDB> <o:complexFile>",
                CITATION_FOLDSEEK_MULTIMER, {
                                           {"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"complexDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"complexFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}
                                   }
        },
        {"createcomplexreport", createmultimerreport, &localPar.createmultimerreport, COMMAND_HIDDEN,
                "", NULL, "", "", CITATION_FOLDSEEK_MULTIMER, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}
        },
        {"expandmultimer", expandmultimer, &localPar.expandmultimer, COMMAND_PREFILTER,
                "Re-prefilter to ensure complete alignment between multimers",
                NULL,
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:alignmentDB> <o:prefilterDB>",
                CITATION_FOLDSEEK_MULTIMER, {
                                        {"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                        {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                        {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                        {"prefilterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::prefilterDb }
                                }
        },
        {"expandcomplex", expandmultimer, &localPar.expandmultimer, COMMAND_HIDDEN,
                "", NULL, "", "", CITATION_FOLDSEEK_MULTIMER, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}
        },
        {"version",              versionstring,        &localPar.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_FOLDSEEK_MULTIMER, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}
};

std::vector<KmerThreshold> externalThreshold = { {Parameters::DBTYPE_AMINO_ACIDS, 7, 197.0, 11.22}};

std::vector<DatabaseDownload> externalDownloads = {
        {
                "Alphafold/UniProt",
                "AlphaFold UniProt Protein Structure Database (including C-alpha, ~700GB download, ~950GB extracted).",
                "Varadi et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Research, (2024)",
                "https://alphafold.ebi.ac.uk/",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "Alphafold/UniProt50-minimal",
                "AlphaFold UniProt Protein Structure Database clustered with MMseqs2 at 50% sequence identity and 90% bidrectional coverage (representative only).",
                "Varadi et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Research, (2024)",
                "https://alphafold.ebi.ac.uk/",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "Alphafold/UniProt50",
                "AlphaFold UniProt Protein Structure Database clustered with MMseqs2 at 50% sequence identity and 90% bidrectional coverage.",
                "Varadi et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Research, (2024)",
                "https://alphafold.ebi.ac.uk/",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "Alphafold/Proteome",
                "AlphaFold Proteomes Protein Structure Database.",
                "Jumper et al. Highly accurate protein structure prediction with AlphaFold. Nature, (2021)",
                "https://alphafold.ebi.ac.uk/",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "Alphafold/Swiss-Prot",
                "AlphaFold Swissprot Protein Structure Database.",
                "Jumper et al. Highly accurate protein structure prediction with AlphaFold. Nature, (2021)",
                "https://alphafold.ebi.ac.uk/",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "ESMAtlas30",
                "ESM Metagenomic Atlas clustered at 30% sequence identity.",
                "Lin et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science, (2023)",
                "https://esmatlas.com",
                false, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "PDB",
                "The Protein Data Bank is the single worldwide archive of structural data of biological macromolecules.",
                "Berman et al. The Protein Data Bank. Nucleic Acids Res, 28(1), 235-242 (2000)",
                "https://www.rcsb.org",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "CATH50",
                "CATH domain database (combined AlphaFold and PDB CATH clustered at 50% seq.id.).",
                "Bordin et al. AlphaFold2 reveals commonalities and novelties in protein structure space for 21 model organisms. Communications Biology, 6, 160 (2023)",
                "https://www.cath.info",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "BFMD",
                "BFMD Big fantastic multimer database (combined multimers from large prediction projects).",
                "Kim W et al. Rapid and Sensitive Protein Complex Alignment with Foldseek-Multimer. bioRxiv (2024)",
                "https://foldseek.steineggerlab.workers.dev/bfmd.version",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "BFVD",
                "BFVD Big fantastic virus database (predicted viral protein structures based on the viral sequence representatives of the UniRef30 clusters",
                "Kim R et al. BFVD - a large repository of predicted viral protein structures. Nucleic Acids Research, gkae1119 (2024)",
                "https://bfvd.steineggerlab.workers.dev",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "ProstT5",
                "Protein language model to predict 3Di directly from sequence.",
                "Heinzinger et al. Bilingual language model for protein sequence and structure. NAR Genomics and Bioinformatics, lqae150 (2024)",
                "https://huggingface.co/Rostlab/ProstT5",
                false, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        }
};