//
// Created by Martin Steinegger on 02.10.18.
//
#include "kmermatcher.h"
#include <mmseqs/src/commons/Timer.h>
#include <mmseqs/src/prefiltering/QueryMatcher.h>
#include <mmseqs/src/prefiltering/ReducedMatrix.h>
#include "StructSubstitutionMatrix.h"


#ifdef OPENMP
#include <omp.h>
#endif

void setStructKmerFilterDefault(Parameters *p) {
    p->spacedKmer = false;
    p->covThr = 0.8;
    p->maskMode = 0;
    p->kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    p->alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
    p->kmersPerSequence = Parameters::CLUST_LINEAR_KMER_PER_SEQ;
}

int structkmermatcher(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    setStructKmerFilterDefault(&par);
    par.parseParameters(argc, argv, command, 2, false, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), par.db1Index.c_str());
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    int querySeqType  =  seqDbr.getDbtype();

    setKmerLengthAndAlphabet(par, seqDbr.getAminoAcidDBSize(), querySeqType);
    std::vector<MMseqsParameter>* params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    Debug(Debug::INFO) << "Database type: " << seqDbr.getDbTypeName() << "\n";

    SubstitutionMatrix structMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    BaseMatrix *subMat = &structMat;
    subMat = new ReducedMatrix(structMat.probMatrix, structMat.subMatrixPseudoCounts, structMat.aa2int, structMat.int2aa, structMat.alphabetSize,  par.alphabetSize, 2.0);


    //seqDbr.readMmapedDataInMemory();
    const size_t KMER_SIZE = par.kmerSize;
    size_t chooseTopKmer = par.kmersPerSequence;

    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = static_cast<size_t>(par.splitMemoryLimit) * 1024;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }
    Debug(Debug::INFO) << "\n";
    size_t totalKmers = computeKmerCount(seqDbr, KMER_SIZE, chooseTopKmer);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter(totalKmers);
    Debug(Debug::INFO) << "Needed memory (" << totalSizeNeeded << " byte) of total memory (" << memoryLimit << " byte)\n";
    // compute splits
    size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
    if (splits > 1) {
        // security buffer
        splits += 1;
    }

    Debug(Debug::INFO) << "Process file into " << splits << " parts\n";
    std::vector<std::string> splitFiles;
    KmerPosition *hashSeqPair = NULL;

    size_t mpiRank = 0;
#ifdef HAVE_MPI
    splits = std::max(static_cast<size_t>(MMseqsMPI::numProc), splits);
    size_t fromSplit = 0;
    size_t splitCount = 1;
    mpiRank = MMseqsMPI::rank;
    // if split size is great than nodes than we have to
    // distribute all splits equally over all nodes
    unsigned int * splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
    memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
    for(size_t i = 0; i < splits; i++){
        splitCntPerProc[i % MMseqsMPI::numProc] += 1;
    }
    for(int i = 0; i < MMseqsMPI::rank; i++){
        fromSplit += splitCntPerProc[i];
    }
    splitCount = splitCntPerProc[MMseqsMPI::rank];
    delete[] splitCntPerProc;

    for(size_t split = fromSplit; split < fromSplit+splitCount; split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation(totalKmers, split, splits, splitFileName, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpiRank == 0){
        for(size_t split = 0; split < splits; split++) {
            std::string splitFileName = par.db2 + "_split_" +SSTR(split);
            splitFiles.push_back(splitFileName);
        }
    }
#else
    for(size_t split = 0; split < splits; split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation(totalKmers, split, splits, splitFileName, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer);
        splitFiles.push_back(splitFileName);
    }
#endif
    if(mpiRank == 0){
        std::vector<char> repSequence(seqDbr.getSize());
        std::fill(repSequence.begin(), repSequence.end(), false);
        // write result
        DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads);
        dbw.open();

        Timer timer;
        if(splits > 1) {
            std::cout << "How many splits: " << splits<<std::endl;
            seqDbr.unmapData();
            mergeKmerFilesAndOutput(seqDbr, dbw, splitFiles, repSequence, par.covMode, par.cov);
        } else {
            writeKmerMatcherResult(seqDbr, dbw, hashSeqPair, totalKmers, repSequence, par.covMode, par.cov, par.threads);
        }
        Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";
        // add missing entries to the result (needed for clustering)

        {
#pragma omp parallel for
            for (size_t id = 0; id < seqDbr.getSize(); id++) {
                char buffer[100];
                int thread_idx = 0;
#ifdef OPENMP
                thread_idx = omp_get_thread_num();
#endif
                if (repSequence[id] == false) {
                    hit_t h;
                    h.pScore = 0;
                    h.diagonal = 0;
                    h.seqId = seqDbr.getDbKey(id);
                    int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                    dbw.writeData(buffer, len, seqDbr.getDbKey(id), thread_idx);
                }
            }
        }
        dbw.close();

    }
    // free memory
    delete subMat;
    if(hashSeqPair){
        delete [] hashSeqPair;
    }
    seqDbr.close();

    return EXIT_SUCCESS;
}

