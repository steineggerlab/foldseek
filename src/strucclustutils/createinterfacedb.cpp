#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "createinterfacedb.sh.h"

int createinterfacedb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string program = par.db2 + ".sh";
    FileUtil::writeFile(program, createinterfacedb_sh, createinterfacedb_sh_len);
    CommandCaller cmd;
    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN", par.filenames.back().c_str());
    cmd.addVariable("CREATESOMEINTERFACEDB_PAR", par.createParameterString(par.createSomeinterfacedb).c_str());
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
