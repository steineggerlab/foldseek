#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "makepaddeddb.sh.h"

int makepaddeddb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string program;
    if (par.db2.find('/') == std::string::npos) {
        program = "./" + par.db2 + ".sh";
    } else {
        program = par.db2 + ".sh";
    }

    FileUtil::writeFile(program, makepaddeddb_sh, makepaddeddb_sh_len);
    CommandCaller cmd;

    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN", par.filenames.back().c_str());

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    if (par.clusterSearch == 1) {
        cmd.addVariable("CLUSEARCH_PAR", "1");
    } else {
        cmd.addVariable("CLUSEARCH_PAR", "0");
    }
    cmd.addVariable("MAKEPADDEDSEQDB_PAR", par.createParameterString(par.makepaddedseqdb).c_str());

    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}