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
    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("LIST", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("CREATESTRUCTSUBDB1_PAR", par.createParameterString(par.createstructsubdb).c_str());
    par.dbIdMode = 0;
    cmd.addVariable("CREATESTRUCTSUBDB2_PAR", par.createParameterString(par.createstructsubdb).c_str());
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}