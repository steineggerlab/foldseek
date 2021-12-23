
#include <vector>
#include <sstream>
#include <string.h>

#include <tmalign/TMalign.h>
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "simd.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <cmath>


//// NN simd version with ungapped alignment to make it faster but sensitiviy suffers :(
//template<const int T>
//simd_int needlemanWunschScoreVec(const simd_int * subQNNi, const simd_int * target_sub_vec,
//                                 const  simd_int *subQ_dist, const  simd_int *subT_dist,
//                                 const simd_int nwGapPenalty, const simd_int vSubMatBias, const simd_int vDistMatBias){
//
//    simd_int prev_sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
//    simd_int sMat_vec[T+1] __attribute__((aligned(ALIGN_INT)));
//    memset(prev_sMat_vec, 0, sizeof(simd_int) * (T+1));
//    simd_int value = simdi_setzero();
//    for(int i = 1; i < T+1; i++){
//        sMat_vec[0]= simdi_setzero();
//
//        // score
//        simd_int scoreLookup = UngappedAlignment::Shuffle(target_sub_vec[i-1], subQNNi[i-1]);
//        scoreLookup = simdi_and(scoreLookup, simdi16_set(0x00FF));
//// do this after the inner loop and use simdi8 TODO
//        scoreLookup = simdi16_sub(scoreLookup, vSubMatBias);
////            for(int a = 0; a < 9; a++){
////                cout << ((short *)&scoreLookup)[a] << " ";
////            }
////            cout << endl;
//        // distance
//        simd_int distLookup = UngappedAlignment::Shuffle(subT_dist[i-1], subQ_dist[i-1]);
//        distLookup = simdi_and(distLookup, simdi16_set(0x00FF));
//        distLookup = simdi16_sub(distLookup, vDistMatBias);
//
//        // add
//        scoreLookup = simdi16_add(scoreLookup, distLookup);
//        value = simdi16_add(value,scoreLookup);
//    }
//
//    return value; /// 4;
//}

int pareunaligner(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> qdbrAA((par.db1).c_str(), (par.db1Index).c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qdbrAA.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> qdbr((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> qcadbr((par.db1+"_ca").c_str(), (par.db1+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qcadbr.open(DBReader<unsigned int>::NOSORT);

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    SubstitutionMatrix subMatAA("/Users/charlotte/CLionProjects/foldseek/lib/mmseqs/data/blosum62.out", 1.4, par.scoreBias);

    DBReader<unsigned int> *tdbr = NULL;
    DBReader<unsigned int> *tdbrAA = NULL;
    DBReader<unsigned int> *tcadbr = NULL;

    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool sameDB = false;

    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
        tdbrAA = &qdbrAA;
        tcadbr = &qcadbr;
    } else {
        tdbr = new DBReader<unsigned int>((par.db2+"_ss").c_str(), (par.db2+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        tdbr->open(DBReader<unsigned int>::NOSORT);

        tdbrAA = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        tdbrAA->open(DBReader<unsigned int>::NOSORT);

        tcadbr = new DBReader<unsigned int>((par.db2+"_ca").c_str(), (par.db2+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        tcadbr->open(DBReader<unsigned int>::NOSORT);
        if (touch) {
            tdbr->readMmapedDataInMemory();
            tcadbr->readMmapedDataInMemory();
        }
    }

    Debug(Debug::INFO) << "Result database: " << par.db3 << "\n";
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug(Debug::INFO) << "Output file: " << par.db4 << "\n";
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    //temporary output file
    Debug::Progress progress(resultReader.getSize());

    // sub. mat needed for query profile
    int8_t * tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat3Di.alphabetSize * 32);
    int8_t * tinySubMatAA = (int8_t*) mem_align(ALIGN_INT, subMatAA.alphabetSize * 32);

    for (int i = 0; i < subMat3Di.alphabetSize; i++) {
        for (int j = 0; j < subMat3Di.alphabetSize; j++) {
            tinySubMat3Di[i * subMat3Di.alphabetSize + j] = subMat3Di.subMatrix[i][j]; // for farrar profile
        }
    }
    for (int i = 0; i < subMatAA.alphabetSize; i++) {
        for (int j = 0; j < subMatAA.alphabetSize; j++) {
            tinySubMatAA[i * subMatAA.alphabetSize + j] = subMatAA.subMatrix[i][j];
        }
    }

    EvalueComputation evaluer(tdbrAA->getAminoAcidDBSize(), &subMatAA);

#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        std::vector<Matcher::result_t> alignmentResult;
        StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, false);

        Sequence qSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);

        char buffer[1024+32768];
        std::string resultBuffer;

        // write output file

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            size_t queryKey = resultReader.getDbKey(id);
            if(*data != '\0') {
                unsigned int queryId = qdbr.getId(queryKey);

                char *querySeqAA = qdbrAA.getData(queryId, thread_idx);
                char *querySeq3Di = qdbr.getData(queryId, thread_idx);

                unsigned int querySeqLen = qdbr.getSeqLen(queryId);
                qSeq3Di.mapSequence(id, queryKey, querySeq3Di, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);
                int32_t maskLen = querySeqLen / 2;

                structureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);

                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                    char * targetSeq3Di = tdbr->getData(targetId, thread_idx);
                    char * targetSeqAA = tdbrAA->getData(targetId, thread_idx);
                    const int targetLen = static_cast<int>(tdbr->getSeqLen(targetId));

                    tSeq3Di.mapSequence(targetId, dbKey, targetSeq3Di, targetLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetLen);
                    if(Util::canBeCovered(par.covThr, par.covMode, qSeq3Di.L, targetLen) == false){
                        continue;
                    }


                    Matcher::result_t res;

                    StructureSmithWaterman::s_align align = structureSmithWaterman.ssw_align(tSeqAA.numSequence, tSeq3Di.numSequence, targetLen, par.gapOpen.values.aminoacid(),
                                                     par.gapExtend.values.aminoacid(), par.alignmentMode,
                                                     par.evalThr, &evaluer, par.covMode, par.covThr, maskLen);

                    unsigned int targetKey = tdbr->getDbKey(targetId);
                    res.score = align.score1;
                    res.eval = align.evalue;

                    res.dbKey = targetKey;
                    res.eval = 0;
                    //Matcher::result_t res = paruenAlign.align(qSeq3Di, tSeq3Di, &subMat3Di, evaluer);

                    if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        alignmentResult.emplace_back(res);
                        passedNum++;
                        rejected = 0;
                    }

                    else {
                        rejected++;
                    }
                }
            }
            if (alignmentResult.size() > 1) {
                SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), Matcher::compareHits);
            }
            for (size_t result = 0; result < alignmentResult.size(); result++) {
                size_t len = Matcher::resultToBuffer(buffer, alignmentResult[result], par.addBacktrace);
                resultBuffer.append(buffer, len);
            }
            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            alignmentResult.clear();
        }
    }

    delete [] tinySubMatAA;
    delete [] tinySubMat3Di;

    dbw.close();
    resultReader.close();
    qdbr.close();
    qcadbr.close();
    if(sameDB == false){
        tdbr->close();
        tcadbr->close();
        delete tdbr;
        delete tcadbr;
    }
    return EXIT_SUCCESS;
}
