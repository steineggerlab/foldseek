#include "Command.h"
#include "CommandDeclarations.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "strucclust";
const char* tool_name = "strucclust";
const char* tool_introduction = "Protein Structure Clustering.";
const char* main_author = "Martin Steinegger (martin.steinegger@mpibpc.mpg.de)";
const char* show_extended_help = NULL;
const char* show_bash_info = NULL;
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
        {"strucclust",             strucclust,            &localPar.strucclust,    COMMAND_MAIN,
                "Assemble protein sequences by iterative greedy overlap assembly.",
                "Extends sequence to the left and right using ungapped alignments.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> ",
                "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
                CITATION_PLASS },
        {"tmalign",      tmalign,      &localPar.tmalign,      COMMAND_HIDDEN,
                "Compute tm-score ",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_MMSEQS2},
        {"convert2statealphabet",    convert2statealphabet,    &localPar.onlythreads,          COMMAND_HIDDEN,
                "Compute state alphabet from dssp input (BLAST3D alphabet)",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:dsspDB> <o:sequenceDB>",
                CITATION_MMSEQS2},
        {"version",              versionstring,        &localPar.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_MMSEQS2}
};
