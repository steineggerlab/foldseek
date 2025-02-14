#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "createinterfacedbs.sh.h"

int createinterfacedbs(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string program = par.db2 + ".sh";
    FileUtil::writeFile(program, createinterfacedbs_sh, createinterfacedbs_sh_len);
    CommandCaller cmd;
    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN", par.filenames.back().c_str());
    //TODO: about parameters
    
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}