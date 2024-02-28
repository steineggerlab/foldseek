#include <cassert>

#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "easycomplexcluster.sh.h"

void setEasyComplexClusterDefaults(Parameters *p) {
    p->removeTmpFiles = true;
}

void setEasyComplexClusterMustPassAlong(Parameters *p) {
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
}

int easycomplexcluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasyComplexClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyComplexClusterMustPassAlong(&par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("INPUT", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("COMPLEXCLUSTER_PAR", par.createParameterString(par.complexclusterworkflow, true).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("RESULT2REPSEQ_PAR", par.createParameterString(par.result2repseq).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    std::string program = tmpDir + "/easycomplexcluster.sh";
    FileUtil::writeFile(program, easycomplexcluster_sh, easycomplexcluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
