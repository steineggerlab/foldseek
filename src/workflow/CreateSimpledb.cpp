#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "createsimpledb.sh.h"


int createsimpledb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    // TODO: hide param dbtype.
    // par.PARAM_DB_TYPE.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    CommandCaller cmd;
    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN", par.filenames.back().c_str());
    cmd.addVariable("CREATESIMPLEDB_PAR", par.createParameterString(par.createsimpledbworkflow).c_str());
    
    std::string program = par.db2 + ".sh";
    FileUtil::writeFile(program, createsimpledb_sh, createsimpledb_sh_len);

    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}