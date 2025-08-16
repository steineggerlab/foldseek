#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "multimercluster.sh.h"

void setMultimerClusterDefaults(LocalParameters *p) {
    p->tmScoreThr = 0.65; // TODO
    p->filtChainTmThr = 0.001; // TODO
    p->filtInterfaceLddtThr = 0.5; // TODO
}
   

void mustsetMultimerCluster(LocalParameters *p) { 
    p->clusteringSetMode = 1;
}

int multimercluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT); //prefilter
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);  //align
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setMultimerClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    mustsetMultimerCluster(&par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    std::cout<<tmpDir.c_str()<<std::endl;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("INPUT", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("MULTIMERSEARCH_PAR", par.createParameterString(par.multimersearchworkflow, true).c_str()); 
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    // cmd.addVariable("VERBCOMPRESS", par.createParameterString(par.verbandcompression).c_str());

    std::string program = tmpDir + "/multimercluster.sh";
    FileUtil::writeFile(program, multimercluster_sh, multimercluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}