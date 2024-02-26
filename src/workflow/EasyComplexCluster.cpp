#include <cassert>

#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "easycomplexcluster.sh.h"

void setEasyComplexClusterDefaults(Parameters *p) {
    //TODO
    // p->PARAM_C = 0.8;
    p->PARAM_COV_MODE = 1;
    p->sensitivity = 4;
    p->PARAM_CLUSTER_MODE = Parameters::GREEDY;
    p->evalThr = 0.001;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->gapOpen = 10;
    p->gapExtend = 1;
}

void setEasyComplexClusterMustPassAlong(Parameters *p) {
    // p->PARAM_C.wasSet = true;
    p->PARAM_E.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_S.wasSet = true;
    par->addBacktrace = true;
    par->PARAM_ADD_BACKTRACE.wasSet = true;

}
int easycomplexcluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);

    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasyComplexSearchDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyComplexClusterMustPassAlong(&par);

    if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM ||
        par.formatAlignmentMode == LocalParameters::FORMAT_ALIGNMENT_PDB_SUPERPOSED  ||
        par.greedyBestHits) {
        needBacktrace = true;
    }
    if (needBacktrace) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed, since they were requested by output format.\n";
        par.addBacktrace = true;
        par.PARAM_ADD_BACKTRACE.wasSet = true;
    }

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

    cmd.addVariable("CLUSTER_MODULE", "complexcluster");
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("COMPLEXSEARCH_PAR", par.createParameterString(par.complexsearchworkflow).c_str());
    cmd.addVariable("COMPLEXCLUSTER_PAR", par.createParameterString(par.complexclusterworkflow).c_str());
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
    return 0;
}
