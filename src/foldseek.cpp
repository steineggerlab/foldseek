#include "DownloadDatabase.h"
#include "Prefiltering.h"
#include "CommandDeclarations.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "foldseek";
const char* tool_name = "foldseek";
const char* tool_introduction = "Foldseek enables fast and sensitive comparisons of large structure sets. It reaches sensitivities similar to state-of-the-art structural aligners while being at least 20,000 times faster.\n\nPlease cite:\nvan Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C.L.M., Söding, J., and Steinegger, M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)";
const char* main_author = "Michel van Kempen, Stephanie Kim, Charlotte Tumescheit, Milot Mirdita, Jeongjae Lee, Cameron L. M. Gilchrist, Johannes Söding, Martin Steinegger";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
const char* index_version_compatible = "fs1";
bool hide_base_commands = true;
bool hide_base_downloads = true;
LocalParameters& localPar = LocalParameters::getLocalInstance();


void updateValdiation(){
    DbValidator::allDb.push_back(LocalParameters::DBTYPE_CA_ALPHA);
    DbValidator::allDb.push_back(LocalParameters::DBTYPE_TMSCORE);
    DbValidator::allDbAndFlat.push_back(LocalParameters::DBTYPE_CA_ALPHA);
    DbValidator::allDbAndFlat.push_back(LocalParameters::DBTYPE_TMSCORE);
}

void (*validatorUpdate)(void) = updateValdiation;


std::vector<struct Command> commands = {
        {"createdb",             structcreatedb,                &localPar.structurecreatedb,    COMMAND_MAIN,
                "Convert PDB/mmCIF/tar[.gz]/DB files to a db",
                "Convert PDB/mmCIF/tar[.gz]/DB files to a db",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]|tar|DB> ... <i:PDB|mmCIF[.gz]|tar|DB> <o:sequenceDB>",
                CITATION_FOLDSEEK, {{"PDB|mmCIF[.gz]|stdin|tar|DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC,
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
        {"compressca",           compressca,             &localPar.onlythreads,           COMMAND_FORMAT_CONVERSION | COMMAND_EXPERT,
                "Create a new compressed C-alpha DB with 16-bit diff encoding where possible from a sequence DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:caDB>",
                CITATION_FOLDSEEK, {{"Db", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::sequenceDb },
                                           {"caDb", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::cadb }}},
        {"convert2pdb",          convert2pdb,             &localPar.onlyverbosity,        COMMAND_FORMAT_CONVERSION,
                "Convert a foldseek structure db to a multi model PDB file",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:Db> <o:pdbFile>",
                CITATION_FOLDSEEK, {{"Db", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"pdbFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"scorecomplex", scorecomplex, &localPar.scorecomplex, COMMAND_ALIGNMENT,
                "Get complex level alignments from alignmentDB",
                "# Get complex level alignments (chain assignments and tm-scores) from alignmentDB.\n"
                "foldseek scorecomplex queryDB targetDB alignmentDB complexDB\n"
                "# simple tsv output format"
                "foldseek createcomplexreport queryDB targetDB complexDB result.tsv"
                "# output files with convertalis"
                "foldseek convertalis queryDB targetDB complexDB result.m8\n\n",
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:queryDb> <i:targetDb> <i:alignmentDB> <o:complexDB>",
                CITATION_FOLDSEEK, {
                                           {"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::NEED_HEADER, &DbValidator::sequenceDb},
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::NEED_HEADER, &DbValidator::sequenceDb},
                                           {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb},
                                           {"complexDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}
                                   }
        },
        {"easy-complexsearch", easycomplexsearch, &localPar.easyscorecomplexworkflow, COMMAND_EASY,
                "Complex level search",
                "# Search a single/multiple PDB file against a set of PDB files and get complex level alignments\n"
                "foldseek easy-complexsearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp\n"
                "# Format output differently\n"
                "foldseek easy-complexsearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp --format-output query,target,qstart,tstart,cigar\n"
                "# Align with TMalign (global)\n"
                "foldseek easy-complexsearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp --alignment-type 1\n"
                "# Skip prefilter and perform an exhaustive alignment (slower but more sensitive)\n"
                "foldseek easy-complexsearch example/1tim.pdb.gz example/8tim.pdb.gz result tmp --exhaustive-search 1\n\n",
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]>|<i:stdin> <i:targetFastaFile[.gz]>|<i:targetDB> <o:outputFileName> <tmpDir>",
                CITATION_FOLDSEEK, {
                                           {"PDB|mmCIF[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &FoldSeekDbValidator::flatfileStdinAndFolder},
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::flatfileAndFolder},
                                           {"outputFileName", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                           {"tempDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}
                                   }
        },
        {"createcomplexreport", createcomplexreport, &localPar.createcomplexreport, COMMAND_FORMAT_CONVERSION,
                "Convert complex DB to tsv format\"",
                "# Create output in tsv format (9 columns):  qComplexName.c_str(), tComplexName.c_str(), qChainString.c_str(), tChainString.c_str(), qTMScore, tTMScore, t, u, assId\n"
                "#  (1,2) identifiers for query and target complex,\n"
                "#  (3,4) chains of query complex and target complex,\n"
                "#  (5,6) tm score based on query and target residue length,\n"
                "#  (8,9) u and t,\n"
                "#  (9) assignment id\n"
                "foldseek convertalis queryDB targetDB complexDB result.tsv\n",
                "Woosub Kim <woosubgo@snu.ac.kr>",
                "<i:queryDb> <i:targetDb> <i:complexDB> <o:complexFile>",
                CITATION_FOLDSEEK, {
                                           {"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                           {"complexDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"complexFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}
                                   }
        },
        {"version",              versionstring,        &localPar.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_FOLDSEEK, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}
};


std::vector<KmerThreshold> externalThreshold = { {Parameters::DBTYPE_AMINO_ACIDS, 7, 197.0, 11.22}};


#include "structdatabases.sh.h"


std::vector<DatabaseDownload> externalDownloads = {
        {
                "Alphafold/UniProt",
                "AlphaFold UniProt Protein Structure Database (including C-alpha, ~700GB download, ~950GB extracted).",
                "Jumper et al. Highly accurate protein structure prediction with AlphaFold. Nature, (2021)",
                "https://alphafold.ebi.ac.uk/",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "Alphafold/UniProt50-minimal",
                "AlphaFold UniProt Protein Structure Database clustered with MMseqs2 at 50% sequence identity and 90% bidrectional coverage (representative only).",
                "Jumper et al. Highly accurate protein structure prediction with AlphaFold. Nature, (2021)",
                "https://alphafold.ebi.ac.uk/",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        },
        {
                "Alphafold/UniProt50",
                "AlphaFold UniProt Protein Structure Database clustered with MMseqs2 at 50% sequence identity and 90% bidrectional coverage.",
                "Jumper et al. Highly accurate protein structure prediction with AlphaFold. Nature, (2021)",
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
                "Lin et al. Evolutionary-scale prediction of atomic level protein structure with a language model. bioRxiv, (2022)",
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
        }
};