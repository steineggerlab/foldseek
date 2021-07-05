#include "CommandDeclarations.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "foldseek";
const char* tool_name = "foldseek";
const char* tool_introduction = "Protein Structure Search and Clustering.";
const char* main_author = "Martin Steinegger (martin.steinegger@snu.ac.kr)";
const char* show_extended_help = NULL;
const char* show_bash_info = NULL;
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();


void updateValdiation(){
    DbValidator::allDb.push_back(LocalParameters::DBTYPE_CA_ALPHA);
    DbValidator::allDb.push_back(LocalParameters::DBTYPE_TMSCORE);
    DbValidator::allDbAndFlat.push_back(LocalParameters::DBTYPE_CA_ALPHA);
    DbValidator::allDbAndFlat.push_back(LocalParameters::DBTYPE_TMSCORE);
}

void (*validatorUpdate)(void) = updateValdiation;


std::vector<struct Command> commands = {
        {"createdb",             createdb,            &localPar.threadsandcompression,    COMMAND_MAIN,
                "Convert PDB/mmCIF files to an db.",
                "Convert PDB/mmCIF files to an db.",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]> <o:sequenceDB>",
                CITATION_MMSEQS2 },
        {"search",               structuresearch,               &localPar.structuresearchworkflow,       COMMAND_MAIN,
                "Sensitive homology search",
                "# Search multiple FASTA against FASTA (like BLASTP, TBLASTN, BLASTX, BLASTN --search-type 3, TBLASTX --search-type 2)\n"
                "mmseqs search queryDB targetDB resultDB tmp\n"
                "mmseqs convertalis queryDB targetDB resultDB result.m8\n\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},


        {"tmalign",      tmalign,      &localPar.tmalign,      COMMAND_HIDDEN,
                "Compute tm-score ",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"tmDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::alignmentDb }}},
        {"aln2tmscore", aln2tmscore,      &localPar.threadsandcompression,      COMMAND_ALIGNMENT,
                "Compute tmscore of an alignment database ",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:alnDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::cadb },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &FoldSeekDbValidator::cadb },
                                   {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                   {"tmDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &FoldSeekDbValidator::tmscore }}},
        {"version",              versionstring,        &localPar.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_MMSEQS2}
};
