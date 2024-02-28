#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "LocalParameters.h"

#include "complexcluster.sh.h"

void setComplexClusterDefaults(Parameters *p) {
    //TODO, parameters for search, filtercomplex, cluster, createresults
    p->covThr = 0.8;
    p->covMode = 1;
    p->clusteringMode = Parameters::GREEDY;
    p->addBacktrace = true;

    // p->sensitivity = 4;
    // p->evalThr = 0.001;
    // p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    // p->gapOpen = 10;
    // p->gapExtend = 1;
}

void setComplexClusterMustPassAlong(Parameters *p) {
    p->PARAM_C.wasSet = true;
    p->PARAM_ADD_BACKTRACE.wasSet = true;
    // p->PARAM_E.wasSet = true;
    // p->PARAM_S.wasSet = true;
    // p->PARAM_ALIGNMENT_MODE.wasSet = true;

}
int complexcluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    // TODO : figure out if commented params needed
    // par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_S.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    // par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setComplexClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setComplexClusterMustPassAlong(&par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;

    cmd.addVariable("COMPLEXSEARCH_PAR", par.createParameterString(par.complexsearchworkflow, true).c_str()); 
    cmd.addVariable("FILTERCOMPLEX_PAR", par.createParameterString(par.filtercomplex).c_str());    
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/complexcluster.sh";
    FileUtil::writeFile(program, complexcluster_sh, complexcluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}