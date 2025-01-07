#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "result2structprofile.sh.h"

int result2structprofile(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string program = par.db4 + ".sh";

    FileUtil::writeFile(program, result2structprofile_sh, result2structprofile_sh_len);

    CommandCaller cmd;
    cmd.addVariable("OUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("RESULT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN2", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("IN1", par.filenames.back().c_str());
    par.filenames.pop_back();

    par.scoringMatrixFile = "blosum62.out";
    par.seedScoringMatrixFile = "blosum62.out";
    cmd.addVariable("PROFILE_PAR", par.createParameterString(par.result2structprofile).c_str());
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
    cmd.addVariable("PROFILE_SS_PAR", par.createParameterString(par.result2structprofile).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
