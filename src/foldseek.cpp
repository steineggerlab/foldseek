#include "Command.h"
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

std::vector<struct Command> commands = {
        {"convert2db",             convert2db,            &localPar.threadsandcompression,    COMMAND_MAIN,
                "Convert PDB/mmCIF files to an db.",
                "Convert PDB/mmCIF files to an db.",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]> <o:sequenceDB>",
                CITATION_MMSEQS2 },
        {"tmalign",      tmalign,      &localPar.tmalign,      COMMAND_HIDDEN,
                "Compute tm-score ",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_MMSEQS2},
        {"convert2statealphabet",    convert2statealphabet,    &localPar.onlythreads,          COMMAND_HIDDEN,
                "Compute state alphabet from dssp input (BLAST3D alphabet)",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:dsspDB> <o:sequenceDB>",
                CITATION_MMSEQS2},
        {"version",              versionstring,        &localPar.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_MMSEQS2}
};
