#include "Parameters.h"
#include "Util.h"
#include "DBWriter.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "structurecluster.sh.h"
#include <cassert>
#include <LocalParameters.h>

void setStructureClusterWorkflowDefaults(LocalParameters *p) {
    p->spacedKmer = true;
    p->covThr = 0.8;
    p->evalThr = 0.01;
    p->sortByStructureBits = 0;
    p->maxResListLen = 1000;
    p->kmersPerSequence = 300;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->maskMode = 0;
    p->compBiasCorrection = 0;
    p->gapOpen = 10;
    p->gapExtend = 1;
}

//TODO this makes no sense for structures
float setAutomaticStructureClusterThreshold(float seqId) {
    float sens;
    if (seqId <= 0.3) {
        sens = 8;
    } else if (seqId > 0.8) {
        sens = 1.0;
    } else {
        sens = 1.0 + (1.0 * (0.7 - seqId) * 10);
    }
    return sens;
}


void setStructuralClusterAutomagicParameters(Parameters& par) {
    if (par.PARAM_NO_COMP_BIAS_CORR.wasSet == false && par.seqIdThr >= 0.7) {
        par.compBiasCorrection = 0;
        par.PARAM_NO_COMP_BIAS_CORR.wasSet = true;
    }

    if (par.PARAM_MIN_DIAG_SCORE.wasSet == false && par.seqIdThr >= 0.7) {
        par.minDiagScoreThr = 60;
        par.PARAM_MIN_DIAG_SCORE.wasSet = true;
    }

    if (par.PARAM_S.wasSet == false) {
        par.sensitivity = setAutomaticStructureClusterThreshold(par.seqIdThr);
        par.PARAM_S.wasSet = true;
        Debug(Debug::INFO) << "Set cluster sensitivity to -s " << par.sensitivity << "\n";
    }

    const bool nonsymetric = (par.covMode == Parameters::COV_MODE_TARGET || par.covMode == Parameters::COV_MODE_QUERY);
    if (par.PARAM_CLUSTER_MODE.wasSet == false) {
        if (nonsymetric) {
            par.clusteringMode = Parameters::GREEDY_MEM;
        } else {
            par.clusteringMode = Parameters::SET_COVER;
        }
        par.PARAM_CLUSTER_MODE.wasSet = true;
        Debug(Debug::INFO) << "Set cluster mode " << ((par.clusteringMode == Parameters::GREEDY_MEM) ? "GREEDY MEM" : "SET COVER") << "\n";
    }
    if (nonsymetric && par.clusteringMode != Parameters::GREEDY && par.clusteringMode != Parameters::GREEDY_MEM) {
        Debug(Debug::WARNING) << "Combining cluster mode " << par.clusteringMode
                              << " in combination with coverage mode " << par.covMode << " can produce wrong results.\n"
                              << "Please use --cov-mode 2\n";
    }
    if (par.singleStepClustering == false && par.clusteringMode == Parameters::CONNECTED_COMPONENT) {
        Debug(Debug::WARNING) << "Connected component clustering produces less clusters in a single step clustering.\n"
                              << "Please use --single-step-cluster";
    }
    if (par.PARAM_CLUSTER_STEPS.wasSet == false) {
        par.clusterSteps = 3;
        par.PARAM_CLUSTER_STEPS.wasSet = true;
        Debug(Debug::INFO) << "Set cluster iterations to " << par.clusterSteps << "\n";
    }
}

int structurecluster(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setStructureClusterWorkflowDefaults(&par);
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_S.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.parseParameters(argc, argv, command, false, 0, 0);
    par.printParameters(command.cmd, argc, argv, *command.params);

    setStructuralClusterAutomagicParameters(par);

    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.clusterworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    std::string alnParam;
    if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_TMALIGN) {
        cmd.addVariable("ALIGNMENT_ALGO", "tmalign");
        alnParam = par.createParameterString(par.tmalign);
    } else if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA || par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI) {
        cmd.addVariable("ALIGNMENT_ALGO", "structurealign");
        alnParam = par.createParameterString(par.structurealign);
    }

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("MERGECLU_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("VERBOSITYANDCOMPRESS", par.createParameterString(par.threadsandcompression).c_str());

    // Linclust parameter
    par.includeIdentity = true;
    cmd.addVariable("STRUCTURERESCOREDIAGONAL_PAR", par.createParameterString(par.structurerescorediagonal).c_str());
    //par.alphabetSize = 14;
    //par.kmerSize = 10;
    //par.spacedKmer = 1;
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());

    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    cmd.addVariable("ALIGNMENT_PAR", alnParam.c_str());
    cmd.addVariable("RUN_LINCLUST", "1");

    if (par.singleStepClustering == false) {
        // save some values to restore them later
        float targetSensitivity = par.sensitivity;
        int kmerSize = par.kmerSize;
        par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
        int maskMode = par.maskMode;
        par.maskMode = 0;
        cmd.addVariable("LINCLUST_PAR", par.createParameterString(par.linclustworkflow).c_str());
        par.kmerSize = kmerSize;
        par.maskMode = maskMode;
        // 1 is lowest sens
        par.sensitivity = ((par.clusterSteps - 1) == 0) ? par.sensitivity : 1;
        int minDiagScoreThr = par.minDiagScoreThr;
        par.minDiagScoreThr = 0;
        par.diagonalScoring = 0;
        par.compBiasCorrection = 0;
        size_t maxResListLen = par.maxResListLen;
        par.maxResListLen = maxResListLen / 10;
        cmd.addVariable("PREFILTER0_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT0_PAR", alnParam.c_str());
        cmd.addVariable("CLUSTER0_PAR", par.createParameterString(par.clust).c_str());
        par.diagonalScoring = 1;
        par.compBiasCorrection = 1;
        par.minDiagScoreThr = minDiagScoreThr;
        float sensStepSize = (targetSensitivity - 1) / (static_cast<float>(par.clusterSteps) - 1);
        for (int step = 1; step < par.clusterSteps; step++) {
            par.sensitivity = 1.0 + sensStepSize * step;
            if(step == (par.clusterSteps - 1)) {
                par.maxResListLen = maxResListLen;
            } else {
                par.maxResListLen = maxResListLen / (10 / (step + 1));
            }

            par.compBiasCorrectionScale = 0.15;
            cmd.addVariable(std::string("PREFILTER" + SSTR(step) + "_PAR").c_str(), par.createParameterString(par.prefilter).c_str());
            par.compBiasCorrectionScale = 0.5;
            cmd.addVariable(std::string("ALIGNMENT" + SSTR(step) + "_PAR").c_str(), alnParam.c_str());
            cmd.addVariable(std::string("CLUSTER" + SSTR(step) + "_PAR").c_str(), par.createParameterString(par.clust).c_str());
        }
        cmd.addVariable("RUN_ITERATIVE", "1");
        cmd.addVariable("STEPS", SSTR(par.clusterSteps).c_str());
        cmd.addVariable("THREADSANDCOMPRESS", par.createParameterString(par.threadsandcompression).c_str());
        cmd.addVariable("VERBCOMPRESS", par.createParameterString(par.verbandcompression).c_str());
        // correct for cascading clustering errors
        if (par.clusterReassignment) {
            cmd.addVariable("REASSIGN", "TRUE");
        }
        int swapedCovMode = Util::swapCoverageMode(par.covMode);
        int tmpCovMode = par.covMode;
        par.covMode = swapedCovMode;
        cmd.addVariable("PREFILTER_REASSIGN_PAR", par.createParameterString(par.prefilter).c_str());
        par.covMode = tmpCovMode;
        cmd.addVariable("ALIGNMENT_REASSIGN_PAR", par.createParameterString(par.structurealign).c_str());
        cmd.addVariable("MERGEDBS_PAR", par.createParameterString(par.mergedbs).c_str());

        std::string program = tmpDir + "/clustering.sh";
        FileUtil::writeFile(program, structurecluster_sh, structurecluster_sh_len);
        cmd.execProgram(program.c_str(), par.filenames);
    } else {
        // same as above, clusthash needs a smaller alphabetsize
        float seqIdThr = par.seqIdThr;
        par.seqIdThr = (float) Parameters::CLUST_HASH_DEFAULT_MIN_SEQ_ID / 100.0f;
        cmd.addVariable("DETECTREDUNDANCY_PAR", par.createParameterString(par.clusthash).c_str());
        par.seqIdThr = seqIdThr;
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        std::string program = tmpDir + "/clustering.sh";
        FileUtil::writeFile(program, structurecluster_sh, structurecluster_sh_len);
        cmd.execProgram(program.c_str(), par.filenames);
    }

    // Unreachable
    assert(false);
    return EXIT_FAILURE;
}

