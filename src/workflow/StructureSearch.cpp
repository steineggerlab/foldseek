#include <cassert>
#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "PrefilteringIndexReader.h"
#include "structuresearch.sh.h"

void setStructureSearchWorkflowDefaults(LocalParameters *p) {
    p->maskMode = 0;
    p->maskProb = 0.99995;
    p->sensitivity = 9.5;
    p->gapOpen = 10;
    p->gapExtend = 1;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->removeTmpFiles = true;
}

void setStructureSearchMustPassAlong(LocalParameters *p) {
    p->PARAM_MASK_RESIDUES.wasSet = true;
    p->PARAM_MASK_PROBABILTY.wasSet = true;
    p->PARAM_NO_COMP_BIAS_CORR.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_GAP_OPEN.wasSet = true;
    p->PARAM_GAP_EXTEND.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
}

int structuresearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    setStructureSearchWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setStructureSearchMustPassAlong(&par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    std::string target = par.filenames.back().c_str();
    cmd.addVariable("TARGET_PREFILTER", (target+"_ss").c_str());
    par.filenames.pop_back();
    std::string query = par.filenames.back().c_str();
    cmd.addVariable("QUERY_PREFILTER", (query+"_ss").c_str());

    const bool isIndex = PrefilteringIndexReader::searchForIndex(target).empty() == false;
    cmd.addVariable("INDEXEXT", isIndex ? ".idx" : NULL);
    par.compBiasCorrectionScale = 0.15;
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
    par.compBiasCorrectionScale = 0.5;
    if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI){
        cmd.addVariable("ALIGNMENT_ALGO", "align");
        cmd.addVariable("QUERY_ALIGNMENT", (query+"_ss").c_str());
        cmd.addVariable("TARGET_ALIGNMENT", (target+"_ss").c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
    }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN){
        cmd.addVariable("ALIGNMENT_ALGO", "tmalign");
        cmd.addVariable("QUERY_ALIGNMENT", query.c_str());
        cmd.addVariable("TARGET_ALIGNMENT", target.c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.tmalign).c_str());
    }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA){
        cmd.addVariable("ALIGNMENT_ALGO", "structurealign");
        cmd.addVariable("QUERY_ALIGNMENT", query.c_str());
        cmd.addVariable("TARGET_ALIGNMENT", target.c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
    }
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    std::string program = tmpDir + "/structuresearch.sh";
    FileUtil::writeFile(program, structuresearch_sh, structuresearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    // Should never get here
    return EXIT_FAILURE;
}
