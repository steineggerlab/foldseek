#include <cassert>

#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "easycomplexcluster.sh.h"

void setEasyComplexClusterDefaults(Parameters *p) {
    //TODO, parameters for search, filtercomplex, cluster, createresults
    p->PARAM_C = 0.8;
    p->PARAM_COV_MODE = 1;
    p->sensitivity = 4;
    p->PARAM_CLUSTER_MODE = Parameters::GREEDY;
    p->evalThr = 0.001;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->gapOpen = 10;
    p->gapExtend = 1;
}

void setEasyComplexClusterMustPassAlong(Parameters *p) {
    p->PARAM_C.wasSet = true;
    p->PARAM_E.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->addBacktrace = true;
    p->PARAM_ADD_BACKTRACE.wasSet = true;

}
int easycomplexcluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);

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

    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("INPUT", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("CLUSTER_MODULE", "filtercomplex");
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("COMPLEXSEARCH_PAR", par.createParameterString(par.complexsearchworkflow).c_str());
    cmd.addVariable("FILTERCOMPLEX_PAR", par.createParameterString(par.filtercomplexworkflow).c_str());
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
