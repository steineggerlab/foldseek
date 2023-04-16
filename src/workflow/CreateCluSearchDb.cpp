#include "Parameters.h"
#include "Util.h"
#include "CommandCaller.h"
#include "FileUtil.h"
#include "createclusearchdb.sh.h"
#include <cassert>
#include <LocalParameters.h>

int createclusearchdb(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.removeTmpFiles = true;
    par.parseParameters(argc, argv, command, false, 0, 0);
    par.printParameters(command.cmd, argc, argv, *command.params);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.clusterworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    std::string alnParam;
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("VERBOSITYANDTHREADS", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("PROFILE_PAR", par.createParameterString(par.result2profile).c_str());
    cmd.addVariable("INPUTDB", par.db1.c_str());
    cmd.addVariable("CLUSTERDB", par.db2.c_str());
    cmd.addVariable("RESULTDB", par.db3.c_str());

    par.pca = 1.4;
    par.pcb = 1.5;
    par.scoringMatrixFile = "3di.out";
    par.seedScoringMatrixFile = "3di.out";
    par.maskProfile = 0;
    par.compBiasCorrection = 0;
    if(par.PARAM_E_PROFILE.wasSet == false){
        par.evalProfile = 0.1;
        par.evalThr = 0.1;
    }
    cmd.addVariable("PROFILE_SS_PAR", par.createParameterString(par.result2profile).c_str());
    std::string program = tmpDir + "/createclusearchdb.sh";
    FileUtil::writeFile(program, createclusearchdb_sh, createclusearchdb_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);
    // Unreachable
    assert(false);
    return EXIT_FAILURE;
}

