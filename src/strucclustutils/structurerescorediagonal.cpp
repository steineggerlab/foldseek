#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"
#include "StructureUtil.h"
#include "DistanceCalculator.h"
#include "QueryMatcher.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "LDDT.h"
#include "structureto12st.h"

#ifdef OPENMP
#include <omp.h>
#endif
#include <memory>



template<typename T>
static DistanceCalculator::LocalAlignment ungappedAlignment(const T *seq3Di1,
                                                            const T *seqAA1,
                                                            const T *seq3Di2,
                                                            const T *seqAA2,
                                                            const unsigned int length,
                                                            short **sub3DiMat,
                                                            short **subAAMat) {
    int maxScore = 0;
    int maxEndPos = 0;
    int maxStartPos = 0;
    int minPos = -1;
//        int maxMinPos = 0;
    int score = 0;
    for(unsigned int pos = 0; pos < length; pos++){
        int curr3Di = sub3DiMat[static_cast<int>(seq3Di1[pos])][static_cast<int>(seq3Di2[pos])];
        int currAAi = subAAMat[static_cast<int>(seqAA1[pos])][static_cast<int>(seqAA2[pos])];
        score = curr3Di + currAAi + score;
        const bool isMinScore = (score <= 0);
        score =  (isMinScore) ? 0 : score;
        minPos = (isMinScore) ? pos : minPos;
        const bool isNewMaxScore = (score > maxScore);
        maxEndPos = (isNewMaxScore) ? pos : maxEndPos;
        maxStartPos = (isNewMaxScore) ? minPos + 1 : maxStartPos;
        maxScore = (isNewMaxScore)? score : maxScore;
    }
    return DistanceCalculator::LocalAlignment(maxStartPos, maxEndPos, maxScore);
}

Matcher::result_t ungappedAlignStructure(Sequence & qSeqAA, Sequence & qSeq3Di, Sequence & qRevSeqAA, Sequence & qRevSeq3Di,
                                         Sequence & tSeqAA, Sequence & tSeq3Di, int diagonal, SubstitutionMatrix & subAAMat,
                                         SubstitutionMatrix & sub3DiMat, EvalueNeuralNet & evaluer,
                                         std::pair<double, double> muLambda, std::string & backtrace, Parameters & par) {
    DistanceCalculator::LocalAlignment res;
    float seqId = 0.0;
    int32_t score = 0;
    backtrace.clear();
    unsigned int minDistToDiagonal = abs(diagonal);
    res.distToDiagonal = minDistToDiagonal;
    res.diagonal = diagonal;
    if (diagonal >= 0 && minDistToDiagonal < static_cast<unsigned int>(qSeqAA.L)) {
        unsigned int minSeqLen = std::min(static_cast<unsigned int>(tSeqAA.L), static_cast<unsigned int>(qSeqAA.L) - minDistToDiagonal);
        res.diagonalLen = minSeqLen;
        DistanceCalculator::LocalAlignment tmp = ungappedAlignment(qSeq3Di.numSequence + minDistToDiagonal,
                                                                   qSeqAA.numSequence + minDistToDiagonal,
                                                                   tSeq3Di.numSequence,
                                                                   tSeqAA.numSequence,
                                                                   minSeqLen, sub3DiMat.subMatrix, subAAMat.subMatrix);

        res.score = tmp.score;
        res.startPos = tmp.startPos;
        res.endPos = tmp.endPos;

        DistanceCalculator::LocalAlignment revTmp = ungappedAlignment(qRevSeq3Di.numSequence + minDistToDiagonal,
                                                                      qRevSeqAA.numSequence + minDistToDiagonal,
                                                                      tSeq3Di.numSequence,
                                                                      tSeqAA.numSequence,
                                                                      minSeqLen, sub3DiMat.subMatrix, subAAMat.subMatrix);

        score = static_cast<int32_t>(tmp.score) - static_cast<int32_t>(revTmp.score);


    } else if (diagonal < 0 && minDistToDiagonal < static_cast<unsigned int>(tSeqAA.L)) {
        unsigned int minSeqLen = std::min(tSeqAA.L - minDistToDiagonal, static_cast<unsigned int>(qSeqAA.L));
        res.diagonalLen = minSeqLen;
        DistanceCalculator::LocalAlignment tmp = ungappedAlignment(qSeq3Di.numSequence, qSeqAA.numSequence,
                                                                   tSeq3Di.numSequence + minDistToDiagonal,
                                                                   tSeqAA.numSequence + minDistToDiagonal,
                                                                   minSeqLen,sub3DiMat.subMatrix, subAAMat.subMatrix);
        res.score = tmp.score;
        res.startPos = tmp.startPos;
        res.endPos = tmp.endPos;

        DistanceCalculator::LocalAlignment revTmp = ungappedAlignment(qRevSeq3Di.numSequence, qSeqAA.numSequence,
                                                                      qRevSeqAA.numSequence + minDistToDiagonal,
                                                                      tSeqAA.numSequence + minDistToDiagonal,
                                                                      minSeqLen, sub3DiMat.subMatrix, subAAMat.subMatrix);


        score = static_cast<int32_t>(tmp.score) - static_cast<int32_t>(revTmp.score);

    }

    double evalue = evaluer.computeEvalueCorr(score, muLambda.first, muLambda.second);

    unsigned int distanceToDiagonal = res.distToDiagonal;
    int diagonalLen = res.diagonalLen;
    float targetCov = static_cast<float>(diagonalLen) / static_cast<float>(tSeqAA.L);
    float queryCov = static_cast<float>(diagonalLen) / static_cast<float>(qSeqAA.L);

    Matcher::result_t result;
    int qStartPos, qEndPos, dbStartPos, dbEndPos;
    // -1 since diagonal is computed from sequence Len which starts by 1
    if (diagonal >= 0) {
        qStartPos = res.startPos + distanceToDiagonal;
        qEndPos = res.endPos + distanceToDiagonal;
        dbStartPos = res.startPos;
        dbEndPos = res.endPos;
    } else {
        qStartPos = res.startPos;
        qEndPos = res.endPos;
        dbStartPos = res.startPos + distanceToDiagonal;
        dbEndPos = res.endPos + distanceToDiagonal;
    }
    unsigned int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, dbStartPos, dbEndPos);
    if (par.addBacktrace) {
        backtrace.append(alnLength, 'M');
    }

    queryCov = SmithWaterman::computeCov(qStartPos, qEndPos, qSeqAA.L);
    targetCov = SmithWaterman::computeCov(dbStartPos, dbEndPos, tSeqAA.L);

    bool hasLowerCoverage = !(Util::hasCoverage(par.covThr, par.covMode, queryCov, targetCov));
    if(hasLowerCoverage){
        return Matcher::result_t(UINT_MAX, score, queryCov, targetCov, seqId, evalue, alnLength,
                                 qStartPos, qEndPos, qSeqAA.L, dbStartPos, dbEndPos, tSeqAA.L, backtrace);
    }
    bool hasLowerEvalue = evalue > par.evalThr;
    if (hasLowerEvalue) {
        return Matcher::result_t(UINT_MAX, score, queryCov, targetCov, seqId, evalue, alnLength,
                                 qStartPos, qEndPos, qSeqAA.L, dbStartPos, dbEndPos, tSeqAA.L, backtrace);
    }else{
        int idCnt = 0;
        for (int i = qStartPos; i <= qEndPos; i++) {
            char qLetter = qSeqAA.numSequence[i];
            char tLetter = tSeqAA.numSequence[dbStartPos + (i - qStartPos)];
            idCnt += (qLetter == tLetter) ? 1 : 0;
        }
        seqId = Util::computeSeqId(par.seqIdMode, idCnt, qSeqAA.L, tSeqAA.L, alnLength);
    }
    return Matcher::result_t(tSeqAA.getDbKey(), score, queryCov, targetCov, seqId, evalue, alnLength,
                             qStartPos, qEndPos, qSeqAA.L, dbStartPos, dbEndPos, tSeqAA.L, backtrace);
}


int structureungappedalign(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbrAA(par.db1, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader qdbr3Di(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);

    IndexReader *t3DiDbr = NULL;
    IndexReader *tAADbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &qdbr3Di;
        tAADbr = &qdbrAA;
    } else {
        tAADbr = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    }

    bool target3Di12St = StructureUtil::is3Di12StDb(t3DiDbr->sequenceReader->getDbtype());
    bool query3Di12St  = StructureUtil::is3Di12StDb(qdbr3Di.sequenceReader->getDbtype());

    bool needTMaligner = (par.tmScoreThr > 0);
    bool needLDDT = (par.lddtThr > 0);

    IndexReader *qcadbr = NULL;
    IndexReader *tcadbr = NULL;
    if(needTMaligner || needLDDT){
        par.addBacktrace = 1;
        qcadbr = new IndexReader(
                par.db1,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca");
        if (sameDB) {
            tcadbr = qcadbr;
        } else {
            tcadbr = new IndexReader(
                    par.db2,
                    par.threads,
                    IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                    touch ? IndexReader::PRELOAD_INDEX : 0,
                    DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                    "_ca"
            );
        }
    }


    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    std::string mat12st;
    if (query3Di12St || target3Di12St) {
        for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
            if (par.substitutionMatrices[i].name == "12st.out") {
                std::string matrixData((const char *)par.substitutionMatrices[i].subMatData,
                                       par.substitutionMatrices[i].subMatDataLen);
                std::string matrixName = par.substitutionMatrices[i].name;
                char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
                mat12st.assign(serializedMatrix);
                free(serializedMatrix);
                break;
            }
        }
        if (mat12st.empty()) {
            Debug(Debug::ERROR) << "Cannot find 12st substitution matrix\n";
            EXIT(EXIT_FAILURE);
        }
    }
    float aaFactor = (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA) ? 1.4 : 0.0;
    SubstitutionMatrix subMatAA(blosum.c_str(), aaFactor, par.scoreBias);
    SubstitutionMatrix *subMat12St = NULL;
    if (query3Di12St || target3Di12St) {
        subMat12St = new SubstitutionMatrix(mat12st.c_str(), par.submat12stScale, par.scoreBias);
    }
    //temporary output file
    Debug::Progress progress(resultReader.getSize());

    // sub. mat needed for query profile
    int8_t * tinySubMatAA = (int8_t*) mem_align(ALIGN_INT, subMatAA.alphabetSize * 32);
    int8_t * tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat3Di.alphabetSize * 32);

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
    int8_t * tinySubMat12St = NULL;
    if (query3Di12St || target3Di12St) {
        // SSW uses subMat3Di.alphabetSize as stride for all matrices including 12st
        int stride = subMat3Di.alphabetSize;
        tinySubMat12St = (int8_t*) mem_align(ALIGN_INT, stride * stride * sizeof(int8_t));
        memset(tinySubMat12St, 0, stride * stride * sizeof(int8_t));
        for (int i = 0; i < subMat12St->alphabetSize; i++) {
            for (int j = 0; j < subMat12St->alphabetSize; j++) {
                tinySubMat12St[i * stride + j] = subMat12St->subMatrix[i][j];
            }
        }
    }

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        EvalueNeuralNet evaluer(tAADbr->sequenceReader->getAminoAcidDBSize(), &subMat3Di);
        std::vector<Matcher::result_t> alignmentResult;
        StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, NULL, NULL);
        StructureSmithWaterman reverseStructureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, NULL, NULL);

        Sequence qSeqAA(par.maxSeqLen, qdbrAA.getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, query3Di12St ? Parameters::DBTYPE_AMINO_ACIDS : qdbr3Di.getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence qRevSeqAA(par.maxSeqLen, qdbrAA.getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qRevSeq3Di(par.maxSeqLen, query3Di12St ? Parameters::DBTYPE_AMINO_ACIDS : qdbr3Di.getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        std::unique_ptr<Sequence> qSeq12St;
        std::unique_ptr<Sequence> tSeq12St;
        std::vector<char> qSeq3Di21Buf;
        std::vector<char> qSeq12StBuf;
        std::vector<char> tSeq3Di21Buf;
        std::vector<char> tSeq12StBuf;
        if (query3Di12St || target3Di12St) {
            if (query3Di12St) {
                qSeq12St.reset(new Sequence(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) subMat12St, 0, false, par.compBiasCorrection));
                qSeq3Di21Buf.reserve(par.maxSeqLen);
                qSeq12StBuf.reserve(par.maxSeqLen);
            }
            if (target3Di12St) {
                tSeq12St.reset(new Sequence(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) subMat12St, 0, false, par.compBiasCorrection));
                tSeq3Di21Buf.reserve(par.maxSeqLen);
                tSeq12StBuf.reserve(par.maxSeqLen);
            }
        }
        TMaligner *tmaligner = NULL;
        if(needTMaligner) {
            tmaligner = new TMaligner(
                    std::max(qdbr3Di.sequenceReader->getMaxSeqLen() + 1, t3DiDbr->sequenceReader->getMaxSeqLen() + 1), false, true, par.exactTMscore);
        }
        LDDTCalculator *lddtcalculator = NULL;
        if(needLDDT) {
            lddtcalculator = new LDDTCalculator(qdbr3Di.sequenceReader->getMaxSeqLen() + 1,  t3DiDbr->sequenceReader->getMaxSeqLen() + 1);
        }

        Coordinate16 qcoords;
        Coordinate16 tcoords;

        std::string backtrace;
        char buffer[1024+32768];
        std::string resultBuffer;

        // write output file

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            size_t queryKey = resultReader.getDbKey(id);
            if(*data != '\0') {
                unsigned int queryId = qdbr3Di.sequenceReader->getId(queryKey);

                char *querySeqAA = qdbrAA.sequenceReader->getData(queryId, thread_idx);
                char *querySeq3Di = qdbr3Di.sequenceReader->getData(queryId, thread_idx);

                unsigned int querySeqLen = qdbr3Di.sequenceReader->getSeqLen(queryId);
                const char *querySeq3Di21 = querySeq3Di;
                if (query3Di12St) {
                    StructureUtil::split3Di12St(querySeq3Di, querySeqLen, qSeq3Di21Buf, qSeq12StBuf, subMat3Di, *subMat12St);
                    querySeq3Di21 = qSeq3Di21Buf.data();
                    qSeq12St->mapSequence(id, queryKey, qSeq12StBuf.data(), querySeqLen);
                }
                qSeq3Di.mapSequence(id, queryKey, querySeq3Di21, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);
                if(needLDDT || needTMaligner){
                    size_t qId = qcadbr->sequenceReader->getId(queryKey);
                    char *qcadata = qcadbr->sequenceReader->getData(qId, thread_idx);
                    size_t qCaLength = qcadbr->sequenceReader->getEntryLen(qId);
                    float* queryCaData = qcoords.read(qcadata, qSeq3Di.L, qCaLength);
                    if(needTMaligner){
                        tmaligner->initQuery(queryCaData, &queryCaData[qSeq3Di.L], &queryCaData[qSeq3Di.L+qSeq3Di.L], NULL, qSeq3Di.L);
                    }
                    if(needLDDT){
                        lddtcalculator->initQuery(qSeq3Di.L, queryCaData, &queryCaData[qSeq3Di.L], &queryCaData[qSeq3Di.L+qSeq3Di.L]);
                    }
                }
                qRevSeq3Di.mapSequence(id, queryKey, querySeq3Di21, querySeqLen);
                qRevSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);
                std::pair<double, double> muLambda = evaluer.predictMuLambda(qSeq3Di.numSequence, qSeq3Di.L);
                structureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA,
                                               query3Di12St ? qSeq12St.get() : NULL, tinySubMat12St);
                qRevSeq3Di.reverse();
                qRevSeqAA.reverse();
                if (query3Di12St) {
                    qSeq12St->reverse();
                }
                reverseStructureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA,
                                                       query3Di12St ? qSeq12St.get() : NULL, tinySubMat12St);
                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    hit_t prefHit = QueryMatcher::parsePrefilterHit(data);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = prefHit.seqId;
                    unsigned int targetId = t3DiDbr->sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                    char * targetSeq3Di = t3DiDbr->sequenceReader->getData(targetId, thread_idx);
                    char * targetSeqAA = tAADbr->sequenceReader->getData(targetId, thread_idx);
                    const int targetSeqLen = static_cast<int>(t3DiDbr->sequenceReader->getSeqLen(targetId));

                    const char *targetSeq3Di21 = targetSeq3Di;
                    if (target3Di12St) {
                        StructureUtil::split3Di12St(targetSeq3Di, targetSeqLen, tSeq3Di21Buf, tSeq12StBuf, subMat3Di, *subMat12St);
                        targetSeq3Di21 = tSeq3Di21Buf.data();
                        tSeq12St->mapSequence(targetId, dbKey, tSeq12StBuf.data(), targetSeqLen);
                    }
                    tSeq3Di.mapSequence(targetId, dbKey, targetSeq3Di21, targetSeqLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);
                    if(Util::canBeCovered(par.covThr, par.covMode, qSeq3Di.L, targetSeqLen) == false){
                        rejected++;
                        continue;
                    }
                    Matcher::result_t res = ungappedAlignStructure(qSeqAA, qSeq3Di, qRevSeqAA, qRevSeq3Di, tSeqAA, tSeq3Di, static_cast<short>(prefHit.diagonal), subMatAA, subMat3Di, evaluer, muLambda, backtrace, par);

                    if(res.dbKey == UINT_MAX){
                        rejected++;
                        continue;
                    }

                    if(needTMaligner) {
                        size_t tId = tcadbr->sequenceReader->getId(res.dbKey);
                        char *tcadata = tcadbr->sequenceReader->getData(tId, thread_idx);
                        size_t tCaLength = tcadbr->sequenceReader->getEntryLen(tId);
                        float* targetCaData = tcoords.read(tcadata, res.dbLen, tCaLength);
                        TMaligner::TMscoreResult tmres = tmaligner->computeTMscore(targetCaData, &targetCaData[res.dbLen], &targetCaData[res.dbLen+res.dbLen], res.dbLen,
                                                                                   res.qStartPos, res.dbStartPos, Matcher::uncompressAlignment(res.backtrace),
                                                                                   TMaligner::normalization(par.tmScoreThrMode, std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos ), res.qLen, res.dbLen));
                        if(tmres.tmscore < par.tmScoreThr){
                            continue;
                        }
                    }
                    if(needLDDT){
                        size_t tId = tcadbr->sequenceReader->getId(res.dbKey);
                        char *tcadata = tcadbr->sequenceReader->getData(tId, thread_idx);
                        size_t tCaLength = tcadbr->sequenceReader->getEntryLen(tId);
                        float* targetCaData = tcoords.read(tcadata, res.dbLen, tCaLength);
                        LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(res.dbLen, res.qStartPos, res.dbStartPos,
                                                                   res.backtrace,
                                                                   targetCaData, &targetCaData[res.dbLen],
                                                                   &targetCaData[res.dbLen+res.dbLen]);
                        if(lddtres.avgLddtScore < par.lddtThr){
                            continue;
                        }
                    }

                    if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        alignmentResult.emplace_back(res);
                        passedNum++;
                        rejected = 0;
                    } else {
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
        if(needTMaligner){
            delete tmaligner;
        }
        if(needLDDT){
            delete lddtcalculator;
        }
    }

    free(tinySubMatAA);
    free(tinySubMat3Di);
    free(tinySubMat12St);

    dbw.close();
    resultReader.close();
    if (sameDB == false) {
        delete t3DiDbr;
        delete tAADbr;
    }
    delete subMat12St;
    return EXIT_SUCCESS;
}
