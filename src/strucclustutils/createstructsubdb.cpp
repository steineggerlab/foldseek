#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "createstructsubdb.sh.h"

int createstructsubdb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string program = par.db3 + ".sh";
    FileUtil::writeFile(program, createstructsubdb_sh, createstructsubdb_sh_len);

    CommandCaller cmd;
    cmd.addVariable("CREATESTRUCTSUBDB_PAR", par.createParameterString(par.createstructsubdb).c_str());
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}