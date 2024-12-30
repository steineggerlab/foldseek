#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "makepaddeddb.sh.h"

// void setmakepaddeddbDefaults(LocalParameters *p) {
//     p->compressed = true;
// }

int makepaddeddb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    std::string program = par.db2 + ".sh";
    FileUtil::writeFile(program, makepaddeddb_sh, makepaddeddb_sh_len);
    
    CommandCaller cmd;
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