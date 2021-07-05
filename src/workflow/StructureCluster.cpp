#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

#include "structurecluster.sh.h"

void setAssemblerWorkflowDefaults(LocalParameters *p) {
    p->spacedKmer = false;
    p->maskMode = 0;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.9;
    p->kmersPerSequence = 60;
    p->alphabetSize = 13;
    p->kmerSize = 14;
}

int strucclust(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setAssemblerWorkflowDefaults(&par);
//    par.overrideParameterDescription((Command &)command, par.PARAM_COV_MODE.uniqid, NULL, NULL, par.PARAM_COV_MODE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_C.uniqid, NULL, NULL, par.PARAM_C.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_MIN_SEQ_ID.uniqid, "Overlap sequence identity threshold [0.0, 1.0]", NULL,  par.PARAM_MIN_SEQ_ID.category);
////    par.overrideParameterDescription((Command &)command, par.PARAM_ORF_MIN_LENGTH.uniqid, "Min codons in orf", "minimum codon number in open reading frames",  par.PARAM_ORF_MIN_LENGTH.category );
//    par.overrideParameterDescription((Command &)command, par.PARAM_NUM_ITERATIONS.uniqid, "Number of assembly iterations [1, inf]", NULL,  par.PARAM_NUM_ITERATIONS.category);
//    par.overrideParameterDescription((Command &)command, par.PARAM_E.uniqid, "Extend sequences if the E-value is below [0.0, inf]", NULL,  par.PARAM_E.category);
//
//    par.overrideParameterDescription((Command &)command, par.PARAM_ID_OFFSET.uniqid, NULL, NULL,  par.PARAM_ID_OFFSET.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_CONTIG_END_MODE.uniqid, NULL, NULL,  par.PARAM_CONTIG_END_MODE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_CONTIG_START_MODE.uniqid, NULL, NULL,  par.PARAM_CONTIG_START_MODE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_ORF_MAX_GAP.uniqid, NULL, NULL,  par.PARAM_ORF_MAX_GAP.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_ORF_START_MODE.uniqid, NULL, NULL,  par.PARAM_ORF_START_MODE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_ORF_FORWARD_FRAMES.uniqid, NULL, NULL,  par.PARAM_ORF_FORWARD_FRAMES.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_ORF_REVERSE_FRAMES.uniqid, NULL, NULL,  par.PARAM_ORF_REVERSE_FRAMES.category | MMseqsParameter::COMMAND_EXPERT);
//
//    par.overrideParameterDescription((Command &)command, par.PARAM_SEQ_ID_MODE.uniqid, NULL, NULL,  par.PARAM_SEQ_ID_MODE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_RESCORE_MODE.uniqid, NULL, NULL,  par.PARAM_RESCORE_MODE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_INCLUDE_ONLY_EXTENDABLE.uniqid, NULL, NULL,  par.PARAM_INCLUDE_ONLY_EXTENDABLE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_KMER_PER_SEQ.uniqid, NULL, NULL,  par.PARAM_KMER_PER_SEQ.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_SORT_RESULTS.uniqid, NULL, NULL,  par.PARAM_SORT_RESULTS.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_TRANSLATION_TABLE.uniqid, NULL, NULL, par.PARAM_TRANSLATION_TABLE.category | MMseqsParameter::COMMAND_EXPERT);
//    par.overrideParameterDescription((Command &)command, par.PARAM_USE_ALL_TABLE_STARTS.uniqid, NULL, NULL, par.PARAM_USE_ALL_TABLE_STARTS.category | MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, 3, true, Parameters::PARSE_VARIADIC);

    CommandCaller cmd;
    std::string tmpPath = par.filenames.back();
    if (FileUtil::directoryExists(tmpPath.c_str()) == false) {
        Debug(Debug::INFO) << "Temporary folder " << tmpPath << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmpPath.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not crate tmp folder " << tmpPath << ".\n";
            return EXIT_FAILURE;
        } else {
            Debug(Debug::INFO) << "Created directory " << tmpPath << "\n";
        }
    }
    size_t hash = par.hashParameter(command.databases, par.filenames, par.strucclust);
    std::string tmpDir = tmpPath + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub folder in temporary directory " << tmpDir << ".\n";
            return EXIT_FAILURE;
        }
    }
    par.filenames.pop_back();
    FileUtil::symlinkAlias(tmpDir, "latest");
    char *p = realpath(tmpDir.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get real path of " << tmpDir << "!\n";
        EXIT(EXIT_FAILURE);
    }
    cmd.addVariable("TMP_PATH", p);
    free(p);

    cmd.addVariable("OUT_FILE", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    // save some values to restore them later
    size_t alphabetSize = par.alphabetSize.aminoacids;
    size_t kmerSize = par.kmerSize;
    bool kmerSizeWasSet = false;
    bool alphabetSizeWasSet = false;
    bool clusterModeSet = false;
    for (size_t i = 0; i < par.strucclust.size(); i++) {
        if (par.strucclust[i]->uniqid == par.PARAM_K.uniqid && par.strucclust[i]->wasSet) {
            kmerSizeWasSet = true;
        }
        if (par.strucclust[i]->uniqid == par.PARAM_ALPH_SIZE.uniqid && par.strucclust[i]->wasSet) {
            alphabetSizeWasSet = true;
        }
    }

    if (kmerSizeWasSet == false) {
        par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    }
    if (alphabetSizeWasSet == false) {
        par.alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
    }

    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
    par.alphabetSize = alphabetSize;
    par.kmerSize = kmerSize;

    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());


    FileUtil::writeFile(tmpDir + "/structurecluster.sh", structurecluster_sh, structurecluster_sh_len);
    std::string program(tmpDir + "/structurecluster.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
