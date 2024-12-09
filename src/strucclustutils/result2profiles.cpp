#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "result2profiles.sh.h"

int result2profiles(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    std::string tmpDir = par.filenames.back(); 
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    std::string program = tmpDir + "/result2profiles.sh";
    FileUtil::writeFile(program, result2profiles_sh, result2profiles_sh_len);
    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("PROFILE_PAR", par.createParameterString(par.result2profiles).c_str());
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
    cmd.addVariable("PROFILE_SS_PAR", par.createParameterString(par.result2profiles).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}