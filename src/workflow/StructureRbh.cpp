#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "LocalParameters.h"
namespace structureRbh{
#include "structurerbh.sh.h"
}

#include <cassert>


void setStructureRbhDefaults(LocalParameters *p) {
    p->sortByStructureBits = 0;
    p->sensitivity = 9.5;
    p->maxResListLen = 1000;
    p->gapOpen = 10;
    p->gapExtend = 1;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->removeTmpFiles = true;
}

int structurerbh(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setStructureRbhDefaults(&par);

    // set a lot of possibly misleading comments to EXPERT mode
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);

    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    // restore threads and verbosity
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);


    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.searchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("SEARCH_A_B_PAR", par.createParameterString(par.structuresearchworkflow).c_str());
    int originalCovMode = par.covMode;
    par.covMode = Util::swapCoverageMode(par.covMode);
    cmd.addVariable("SEARCH_B_A_PAR", par.createParameterString(par.structuresearchworkflow).c_str());
    par.covMode = originalCovMode;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN){
        cmd.addVariable("ALIGNMENT_ALGO", "tmalign");
        par.tmScoreThr = 0.0f;
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.tmalign).c_str());
    }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA){
        cmd.addVariable("ALIGNMENT_ALGO", "structurealign");
        par.evalThr = 100000000;
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
    }

    cmd.addVariable("VERB_COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    std::string program = tmpDir + "/rbh.sh";
    FileUtil::writeFile(program, structureRbh::structurerbh_sh, structureRbh::structurerbh_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
