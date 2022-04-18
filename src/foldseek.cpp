#include "DownloadDatabase.h"
#include "CommandDeclarations.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "foldseek";
const char* tool_name = "foldseek";
const char* tool_introduction = "Foldseek enables fast and sensitive comparisons of large structure sets. It reaches sensitivities similar to state-of-the-art structural aligners while being at least 20,000 times faster.\n\nPlease cite: van Kempen M, Kim S,Tumescheit C, Mirdita M, Söding J, and Steinegger M. Foldseek: fast and accurate protein structure search. bioRxiv, doi:10.1101/2022.02.07.479398 (2021)";
const char* main_author = "Michel van Kempen, Stephanie Kim, Charlotte Tumescheit, Milot Mirdita, Johannes Söding, Martin Steinegger";
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
        {"createdb",             createdb,            &localPar.structurecreatedb,    COMMAND_MAIN,
                "Convert PDB/mmCIF files to an db.",
                "Convert PDB/mmCIF files to an db.",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]> <o:sequenceDB>",
                CITATION_FOLDSEEK, {{"PDB|mmCIF[.gz]|stdin", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfileStdinAndGeneric },
                                          {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"easy-search",          easystructuresearch,           &localPar.easystructuresearchworkflow,   COMMAND_EASY,
                "Sensitive homology search",
                "# Search a single/multiple PDB file against a set of PDB files\n"
                "foldseek easy-search example/d1asha_ example/ result.m8 tmp\n"
                "# Format output differently\n"
                "foldseek easy-search example/d1asha_ example/ result.m8 tmp --format-output query,target,qstart,tstart,cigar\n"
                "# Align with TMalign (global)\n"
                "foldseek easy-search example/d1asha_ example/ result.m8 tmp --alignment-type 1\n\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]>|<i:stdin> <i:targetFastaFile[.gz]>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_FOLDSEEK, {{"fastaFile[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &FoldSeekDbValidator::flatfileStdinAndFolder },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::flatfileAndFolder },
                                          {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-msa",          easymsa,           &localPar.easymsaworkflow,   COMMAND_EASY,
                "Sensitive homology search and build MSAs",
                "# Align a set of PDB files and create a MSA\n"
                "foldseek easy-msa example/d1asha_ result.m8 tmp\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]>|<i:stdin> <o:alignmentFile> <tmpDir>",
                CITATION_FOLDSEEK, {{"PDB|mmCIF[.gz]|stdin", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfileStdinAndGeneric },
                                           {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"search",               structuresearch,               &localPar.structuresearchworkflow,       COMMAND_MAIN,
                "Sensitive homology search",
                "# Search multiple structures (cif,PDB) against structures.\n"
                "foldseek search queryDB targetDB resultDB tmp\n"
                "foldseek convertalis queryDB targetDB resultDB result.m8\n\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
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
        {"structurealign",      structurealign,      &localPar.align,      COMMAND_ALIGNMENT,
                "Compute structural alignment using 3Di alphabet, amino acids and neighborhood information",
                NULL,
                "Charlotte Tumescheit <ch.tumescheit@gmail.com> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"alnDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::alignmentDb }}},
        {"generatetree",      generatetree,      &localPar.align,      COMMAND_ALIGNMENT,
                "Compute guide tree based on Foldseek scores",
                NULL,
                "Cameron Gilchrist <gamcil@snu.ac.kr> & Charlotte Tumescheit <ch.tumescheit@gmail.com> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:tree>",
                CITATION_FOLDSEEK, {{"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                           {"tree", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"traversetree",      traversetree,      &localPar.align,      COMMAND_ALIGNMENT,
                "Traverse guide tree and build progressive MSA",
                NULL,
                "Cameron Gilchrist <gamcil@snu.ac.kr> & Charlotte Tumescheit <ch.tumescheit@gmail.com> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <i:tree> <o:alnDB>",
                CITATION_FOLDSEEK, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"tree", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                           {"alnDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::alignmentDb }}},
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

        {"version",              versionstring,        &localPar.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_FOLDSEEK, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}
};

#include "structdatabases.sh.h"

std::vector<DatabaseDownload> externalDownloads = {
        {
                "Alphafold/Proteome",
                "AlphaFold Protein Structure Database.",
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
                "PDB",
                "The Protein Data Bank is the single worldwide archive of structural data of biological macromolecules.",
                "Berman et al. The Protein Data Bank. Nucleic Acids Res, 28(1), 235-242 (2000)",
                "https://www.rcsb.org",
                true, Parameters::DBTYPE_AMINO_ACIDS, structdatabases_sh, structdatabases_sh_len,
                {}
        }
};
