#include <cassert>

#include "LocalParameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "easymultimersearch.sh.h"

int easymultimersearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.structurecreatedb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }

    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, false, Parameters::PARSE_VARIADIC, 0);
    if(par.PARAM_FORMAT_OUTPUT.wasSet == false){
        par.outfmt = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid";
    }
    par.addBacktrace = true;
    par.PARAM_ADD_BACKTRACE.wasSet = true;
    par.printParameters(command.cmd, argc, argv, *command.params);

    bool needBacktrace = false;
    bool needTaxonomy = false;
    bool needTaxonomyMapping = false;
    bool needLookup = false;
    {
        bool needSequenceDB = false;
        bool need3DiDB = false;
        bool needFullHeaders = false;
        bool needSource = false;
        bool needQCA = false;
        bool needTCA = false;
        bool needTMalign = false;
        bool needLDDT = false;
        LocalParameters::getOutputFormat(
            par.formatAlignmentMode, par.outfmt, needSequenceDB, need3DiDB, needBacktrace, needFullHeaders,
            needLookup, needSource, needTaxonomyMapping, needTaxonomy, needQCA, needTCA, needTMalign, needLDDT
        );
    }

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
    if (needLookup) {
        par.writeLookup = true;
    }

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    CommandCaller cmd;
    if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN){
        cmd.addVariable("COMPLEX_ALIGNMENT_ALGO", "tmalign");
        cmd.addVariable("COMPLEX_ALIGN_PAR", par.createParameterString(par.tmalign).c_str());
    }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA || par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI){
        cmd.addVariable("COMPLEX_ALIGNMENT_ALGO", "structurealign");
        cmd.addVariable("COMPLEX_ALIGN_PAR", par.createParameterString(par.structurealign).c_str());
    }

//    if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN){
//        cmd.addVariable("ALIGNMENT_ALGO", "tmalign");
//        cmd.addVariable("QUERY_ALIGNMENT", query.c_str());
//        cmd.addVariable("TARGET_ALIGNMENT", target.c_str());
//        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.tmalign).c_str());
//        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
//        par.sortByStructureBits = 0;
//        //par.evalThr = 10; we want users to adjust this one. Our default is 10 anyhow.
//        cmd.addVariable("STRUCTUREALIGN_PAR", par.createParameterString(par.structurealign).c_str());
//    }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA || par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI){
//        cmd.addVariable("ALIGNMENT_ALGO", "structurealign");
//        cmd.addVariable("QUERY_ALIGNMENT", query.c_str());
//        cmd.addVariable("TARGET_ALIGNMENT", target.c_str());
//        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.structurealign).c_str());
//    }

    switch(par.prefMode){
        case LocalParameters::PREF_MODE_KMER:
            cmd.addVariable("PREFMODE", "KMER");
            break;
        case LocalParameters::PREF_MODE_UNGAPPED:
            cmd.addVariable("PREFMODE", "UNGAPPED");
            break;
        case LocalParameters::PREF_MODE_EXHAUSTIVE:
            cmd.addVariable("PREFMODE", "EXHAUSTIVE");
            break;
    }
    if(par.exhaustiveSearch){
        cmd.addVariable("PREFMODE", "EXHAUSTIVE");
    }
    cmd.addVariable("NO_REPORT", par.multimerReportMode == 0 ? "TRUE" : NULL);
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("OUTPUT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("TARGET", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("QUERY", par.filenames.back().c_str());
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);
    cmd.addVariable("GPU", par.gpu ? "TRUE" : NULL);
    cmd.addVariable("MAKEPADDEDSEQDB_PAR", par.createParameterString(par.makepaddeddb).c_str());
    par.filenames.pop_back();
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.structurecreatedb).c_str());
    cmd.addVariable("MULTIMERSEARCH_PAR", par.createParameterString(par.multimersearchworkflow, true).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("REPORT_PAR", par.createParameterString(par.createmultimerreport).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    std::string program = tmpDir + "/easymultimersearch.sh";
    FileUtil::writeFile(program, easymultimersearch_sh, easymultimersearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);
    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}