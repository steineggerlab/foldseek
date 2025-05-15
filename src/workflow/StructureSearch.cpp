#include <cassert>
#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "PrefilteringIndexReader.h"
#include "structuresearch.sh.h"
#include "structureiterativesearch.sh.h"
#include "structureprofile.sh.h"


void setStructureSearchWorkflowDefaults(LocalParameters *p) {
    p->kmerSize = 0;
    p->sensitivity = 9.5;
    p->maxResListLen = 1000;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->removeTmpFiles = true;
}

void setStructureSearchMustPassAlong(LocalParameters *p) {
    p->PARAM_K.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_MAX_SEQS.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
}

int structuresearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    setStructureSearchWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setStructureSearchMustPassAlong(&par);
    if((par.alignmentMode == 1 || par.alignmentMode == 2) && par.sortByStructureBits){
        Debug(Debug::WARNING) << "Cannot use --sort-by-structure-bits 1 with --alignment-mode 1 or 2\n";
        Debug(Debug::WARNING) << "Disabling --sort-by-structure-bits\n";
        par.sortByStructureBits = false;
    }


    {
        bool needBacktrace = false;
        bool needTaxonomy = false;
        bool needTaxonomyMapping = false;
        bool needLookup = false;
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

        // check if databases have Calpha coordinates
        if (needQCA) {
            std::string caDB = par.db1 + "_ca.dbtype";
            if (!FileUtil::fileExists(caDB.c_str())) {
                Debug(Debug::ERROR)
                        << "Query database does not contain C-alpha coordinates. Please recreate the database.";
                EXIT(EXIT_FAILURE);
            }
        }
        if (needTCA) {
            std::string caDB = par.db2 + "_ca.dbtype";
            if (!FileUtil::fileExists(caDB.c_str())) {
                Debug(Debug::ERROR)
                        << "Target database does not contain C-alpha coordinates. Please recreate the database.";
                EXIT(EXIT_FAILURE);
            }
        }
    }


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
    double prevEvalueThr = par.evalThr;
    par.evalThr = std::numeric_limits<double>::max();
    cmd.addVariable("UNGAPPEDPREFILTER_PAR", par.createParameterString(par.ungappedprefilter).c_str());
    par.evalThr = prevEvalueThr;
    par.compBiasCorrectionScale = 0.5;

    // GPU can only use the ungapped prefilter
    if (par.gpu == 1 && par.PARAM_PREF_MODE.wasSet == false) {
        par.prefMode = Parameters::PREF_MODE_UNGAPPED;
    }
    
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
   if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN){
        cmd.addVariable("ALIGNMENT_ALGO", "tmalign");
        cmd.addVariable("QUERY_ALIGNMENT", query.c_str());
        cmd.addVariable("TARGET_ALIGNMENT", target.c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.tmalign).c_str());
        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
        par.sortByStructureBits = 0;
        par.tmScoreThrMode = 0.0;
        par.lddtThr = 0.0;
       //par.evalThr = 10; we want users to adjust this one. Our default is 10 anyhow.
        cmd.addVariable("STRUCTUREALIGN_PAR", par.createParameterString(par.structurealign).c_str());
    }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA || par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI){
        cmd.addVariable("ALIGNMENT_ALGO", "structurealign");
        cmd.addVariable("QUERY_ALIGNMENT", query.c_str());
        cmd.addVariable("TARGET_ALIGNMENT", target.c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.structurealign).c_str());
    }
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    std::string program;
    if(par.numIterations > 1 || par.numIterations == 0){
	//par.evalProfile = 0.1;
        double originalEval = par.evalThr;

        par.evalThr = (par.evalThr < par.evalProfile) ? par.evalThr  : par.evalProfile;
        if (par.numIterations == 0) {
            program = tmpDir + "/structureprofile.sh";
            FileUtil::writeFile(program, structureprofile_sh, structureprofile_sh_len);
            par.numIterations = 2;
        } else {
            if(par.PARAM_E_PROFILE.wasSet==false) {
                par.evalThr = 0.001;
            }
            program = tmpDir + "/structureiterativesearch.sh";
            FileUtil::writeFile(program, structureiterativesearch_sh, structureiterativesearch_sh_len);
        }
        for (int i = 0; i < par.numIterations; i++) {
            if (i == (par.numIterations - 1)) {
                par.evalThr = originalEval;
            }
            par.addBacktrace = true;
            par.compBiasCorrectionScale = 0.15;
            cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(),
                            par.createParameterString(par.prefilter).c_str());
            cmd.addVariable(std::string("UNGAPPEDPREFILTER_PAR_" + SSTR(i)).c_str(),
                            par.createParameterString(par.ungappedprefilter).c_str());
            par.compBiasCorrectionScale = 0.5;
            if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI){
                cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.align).c_str());
            }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN){
                cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.tmalign).c_str());
            }else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA){
                cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.structurealign).c_str());
            }
        }
        if(par.clusterSearch == 1){
            cmd.addVariable("MERGERESULTBYSET_PAR", par.createParameterString(par.mergeresultsbyset).c_str());
            cmd.addVariable("EXPAND", "1");
        }

        cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
        //do not used PROFILE
        //cmd.addVariable("PROFILE_PAR", par.createParameterString(par.result2structprofile).c_str());
        cmd.addVariable("VERBOSITY_THREADS_PAR", par.createParameterString(par.threadsandcompression).c_str());
        if(par.PARAM_E_PROFILE.wasSet){
            std::vector<MMseqsParameter*> tmpVec;
            tmpVec.push_back(&par.PARAM_E_PROFILE);
            cmd.addVariable("PROFILE_EVAL",par.createParameterString(tmpVec).c_str() );
        }
        cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());
        cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    }else{
        program = tmpDir + "/structuresearch.sh";
        FileUtil::writeFile(program, structuresearch_sh, structuresearch_sh_len);
    }

    if(par.clusterSearch == 1) {
        if (isIndex){
            // check if we have  SRC_SEQUENCES, SEQUENCES, ALIGNMENT
            DBReader<unsigned int> indexReader( (par.db2+".idx").c_str(),
                                                (par.db2+".idx.index").c_str(),
                                                1, DBReader<unsigned int>::USE_INDEX);
            indexReader.open(DBReader<unsigned int>::NOSORT);
            size_t alignmentIdx = indexReader.getId(PrefilteringIndexReader::ALNINDEX);
            if(alignmentIdx == UINT_MAX){
                Debug(Debug::ERROR)
                        << "Require idx with alignments/cluster for cluster search.";
                EXIT(EXIT_FAILURE);
            }
            indexReader.close();
        }else {
            std::vector<std::string> dbsToCheck = {"_seq", "_seq_ss", "_seq_h"};
            for (size_t i = 0; i < dbsToCheck.size(); i++) {
                std::string db = par.db2 + dbsToCheck[i] + ".dbtype";
                if (!FileUtil::fileExists(db.c_str())) {
                    Debug(Debug::ERROR)
                            << "Require " << db << " database for cluster search.";
                    EXIT(EXIT_FAILURE);
                }
            }
            if (!FileUtil::fileExists((par.db2 + "_clu.dbtype").c_str()) &&
                !FileUtil::fileExists((par.db2 + "_aln.dbtype").c_str())) {
                Debug(Debug::ERROR)
                        << "Require " << par.db2 + "_clu  or " << par.db2 + "_aln database for cluster search.";
                EXIT(EXIT_FAILURE);
            }
        }
        cmd.addVariable("MERGERESULTBYSET_PAR", par.createParameterString(par.mergeresultsbyset).c_str());
        cmd.addVariable("EXPAND", "1");
    }
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    // Should never get here
    return EXIT_FAILURE;
}
