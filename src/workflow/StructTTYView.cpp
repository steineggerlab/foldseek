#include "structty.h"
#include "LocalParameters.h"

int structtyview(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    structty::RunOptions opts;
    opts.input_files.push_back(par.filenames[0]);   // query PDB/CIF
    opts.foldseek_file = par.filenames[1];           // result .m8
    if (par.filenames.size() > 2) {
        opts.foldseek_db = par.filenames[2];         // target DB (optional)
    }

    structty::run(opts);
    return EXIT_SUCCESS;
}
