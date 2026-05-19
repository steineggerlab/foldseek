#include "Parameters.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "structureclusterupdate.sh.h"
#include <cassert>
#include "LocalParameters.h"

void setStructureClusterUpdateDefaults(LocalParameters *p) {
    p->covThr = 0.8;
    p->evalThr = 0.001;
    p->sortByStructureBits = 0;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->compBiasCorrection = 0;
}

void setStructureClusterUpdateMustPassAlong(LocalParameters *p) {
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
}

int structureclusterupdate(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setStructureClusterUpdateDefaults(&par);
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ALT_ALIGNMENT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ_SCALE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_START_SENS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SENS_STEPS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_CLUSTER_REASSIGN.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_NUM_ITERATIONS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);
    setStructureClusterUpdateMustPassAlong(&par);

    std::string tmpDir = par.db6;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.clusterUpdate));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RECOVER_DELETED", par.recoverDeleted ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    int oldVerbosity = par.verbosity;
    par.verbosity = std::min(par.verbosity, 1);
    cmd.addVariable("NOWARNINGS_PAR", par.createParameterString(par.onlyverbosity).c_str());
    par.verbosity = oldVerbosity;

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("DIFF_PAR", par.createParameterString(par.diff).c_str());
    cmd.addVariable("RESULT2REPSEQ_PAR", par.createParameterString(par.result2repseq).c_str());
    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusterworkflow, true).c_str());

    // structural alignment algorithm and parameters
    std::string alnParam;
    if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN) {
        cmd.addVariable("ALIGNMENT_ALGO", "tmalign");
        alnParam = par.createParameterString(par.tmalign);
    } else if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA ||
               par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI) {
        cmd.addVariable("ALIGNMENT_ALGO", "structurealign");
        alnParam = par.createParameterString(par.structurealign);
    } else {
        Debug(Debug::ERROR) << "Unknown alignment type: " << par.alignmentType << "\n";
        EXIT(EXIT_FAILURE);
    }
    cmd.addVariable("ALIGNMENT_PAR", alnParam.c_str());

    // prefilter uses 3Di k-mer lookup
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());

    std::string program = tmpDir + "/structureclusterupdate.sh";
    FileUtil::writeFile(program, structureclusterupdate_sh, structureclusterupdate_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Unreachable
    assert(false);
    return EXIT_FAILURE;
}
