#include <cassert>

#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"

#include "easymultimercluster.sh.h"

void setEasyMultimerClusterDefaults(Parameters *p) {
    //TODO
    p->removeTmpFiles = true;
    p->writeLookup = true;
}

void setEasyMultimerClusterMustPassAlong(Parameters *p) {
    //TODO
    p->clusteringSetMode = 1;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
    p->PARAM_WRITE_LOOKUP.wasSet = true;
}

int easymultimercluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    //TODO
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT); //prefilter
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);  //align
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    par.PARAM_S.addCategory(MMseqsParameter::COMMAND_EXPERT);
    
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasyMultimerClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyMultimerClusterMustPassAlong(&par);

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
    cmd.addVariable("GPU", par.gpu ? "TRUE" : NULL);
    cmd.addVariable("MAKEPADDEDSEQDB_PAR", par.createParameterString(par.makepaddeddb).c_str());
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("MULTIMERCLUSTER_PAR", par.createParameterString(par.multimerclusterworkflow,true).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("CREATESUBDB_PAR", par.createParameterString(par.createsubdb).c_str());
    cmd.addVariable("RESULT2REPSEQ_PAR", par.createParameterString(par.result2repseq).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    std::string program = tmpDir + "/easymultimercluster.sh";
    FileUtil::writeFile(program, easymultimercluster_sh, easymultimercluster_sh_len);

    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
