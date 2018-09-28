#include "LocalParameters.h"
#include "DistanceCalculator.h"
#include "Matcher.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"

#include <limits>
#include <cstdint>
#include <queue>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

class CompareResultBySeqId {
public:
    bool operator() (const Matcher::result_t & r1,const Matcher::result_t & r2) {
//        if(r1.seqId < r2.seqId )
//            return true;
//        if(r2.seqId < r1.seqId )
//            return false;
        if(r1.score < r2.score )
            return true;
        if(r2.score < r1.score )
            return false;
        if(r1.dbKey > r2.dbKey )
            return true;
        if(r2.dbKey > r1.dbKey )
            return false;
        /*  int seqLen1 = r1.qEndPos - r1.qStartPos;
          int seqLen2 = r2.qEndPos - r2.qStartPos;
          if(seqLen1 < seqLen2)
              return true;
          if(seqLen2 < seqLen1 )
              return false;*/
        return false;
    }
};

typedef std::priority_queue<Matcher::result_t, std::vector<Matcher::result_t> , CompareResultBySeqId> QueueBySeqId;
Matcher::result_t selectFragmentToExtend(QueueBySeqId &alignments,
                                             unsigned int queryKey) {
    // results are ordered by score
    while (alignments.empty() == false){
        Matcher::result_t res = alignments.top();
        alignments.pop();
        size_t dbKey = res.dbKey;
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 && res.qStartPos == 0);
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != res.dbLen-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != res.qLen-1);
        const bool isNotIdentity = (dbKey != queryKey);
        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){
            return res;
        }
    }
    return Matcher::result_t(UINT_MAX,0,0,0,0,0,0,0,0,0,0,0,0,"");
}


int doassembly(LocalParameters &par) {
    DBReader<unsigned int> *sequenceDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
    sequenceDbr->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> * alnReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
    alnReader->open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    resultWriter.open();
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0f);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(subMat);

    unsigned char * wasExtended = new unsigned char[sequenceDbr->getSize()];
    std::fill(wasExtended, wasExtended+sequenceDbr->getSize(), 0);

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        std::vector<Matcher::result_t> alignments;
        alignments.reserve(300);
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
            Debug::printProgress(id);
            unsigned int queryId = sequenceDbr->getDbKey(id);
            char *querySeq = sequenceDbr->getData(id);
            unsigned int querySeqLen = sequenceDbr->getSeqLens(id) - 2;
            unsigned int leftQueryOffset = 0;
            unsigned int rightQueryOffset = 0;
            std::string query(querySeq, querySeqLen); // no /n/0
            char *alnData = alnReader->getDataByDBKey(queryId);
            alignments.clear();
            Matcher::readAlignmentResults(alignments, alnData);
            QueueBySeqId alnQueue;
            bool queryCouldBeExtended = false;
            while(alignments.size() > 1){
                bool queryCouldBeExtendedLeft = false;
                bool queryCouldBeExtendedRight = false;
                for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {
                    float scorePerCol = static_cast<float>(alignments[alnIdx].score) / static_cast<float>(alignments[alnIdx].alnLength);
                    float alnLen = static_cast<float>(alignments[alnIdx].alnLength);
                    float ids = static_cast<float>(alignments[alnIdx].seqId) * alnLen;
                    alignments[alnIdx].seqId = ids / (alnLen + 0.5);
                    alignments[alnIdx].score = static_cast<int>(scorePerCol*100);

                    alnQueue.push(alignments[alnIdx]);
                    if (alignments.size() > 1)
                        __sync_or_and_fetch(&wasExtended[sequenceDbr->getId(alignments[alnIdx].dbKey)],
                                            static_cast<unsigned char>(0x40));
                }
                std::vector<Matcher::result_t> tmpAlignments;
                Matcher::result_t besttHitToExtend;
                while ((besttHitToExtend = selectFragmentToExtend(alnQueue, queryId)).dbKey != UINT_MAX) {
                    querySeqLen = query.size();
                    querySeq = (char *) query.c_str();

//                querySeq.mapSequence(id, queryKey, query.c_str());
                    unsigned int targetId = sequenceDbr->getId(besttHitToExtend.dbKey);
                    if (targetId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Could not find targetId  " << besttHitToExtend.dbKey
                                            << " in database " << sequenceDbr->getDataFileName() << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    char *targetSeq = sequenceDbr->getData(targetId);
                    unsigned int targetSeqLen = sequenceDbr->getSeqLens(targetId) - 2;
                    // check if alignment still make sense (can extend the query)
                    if (besttHitToExtend.dbStartPos == 0) {
                        if ((targetSeqLen - (besttHitToExtend.dbEndPos + 1)) <= rightQueryOffset) {
                            continue;
                        }
                    } else if (besttHitToExtend.qStartPos == 0) {
                        if (besttHitToExtend.dbStartPos <= leftQueryOffset) {
                            continue;
                        }
                    }
                    __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x10));
                    int qStartPos, qEndPos, dbStartPos, dbEndPos, score;
                    int diagonal = (leftQueryOffset + besttHitToExtend.qStartPos) - besttHitToExtend.dbStartPos;
                    int dist = std::max(abs(diagonal), 0);
                    if (diagonal >= 0) {
//                    targetSeq.mapSequence(targetId, besttHitToExtend.dbKey, dbSeq);
                        size_t diagonalLen = std::min(targetSeqLen, querySeqLen - abs(diagonal));
                        DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                                querySeq + abs(diagonal),
                                targetSeq, diagonalLen, fastMatrix.matrix);
                        qStartPos = alignment.startPos + dist;
                        qEndPos = alignment.endPos + dist;
                        dbStartPos = alignment.startPos;
                        dbEndPos = alignment.endPos;
                        score = alignment.score;
                    } else {
                        size_t diagonalLen = std::min(targetSeqLen - abs(diagonal), querySeqLen);
                        DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                                querySeq,
                                targetSeq + abs(diagonal),
                                diagonalLen, fastMatrix.matrix);
                        qStartPos = alignment.startPos;
                        qEndPos = alignment.endPos;
                        dbStartPos = alignment.startPos + dist;
                        dbEndPos = alignment.endPos + dist;
                        score = alignment.score;
                    }

                    if (dbStartPos == 0 && qEndPos == (querySeqLen - 1) ) {
                        if(queryCouldBeExtendedRight == true) {
                            float alnLen = qEndPos - qStartPos;
                            float scorePerCol = static_cast<float>(score) / (alnLen+0.5);
                            besttHitToExtend.score = static_cast<int>(scorePerCol*100);
                            tmpAlignments.push_back(besttHitToExtend);
                            continue;
                        }
                        size_t dbFragLen = (targetSeqLen - dbEndPos) - 1; // -1 get not aligned element
                        std::string fragment = std::string(targetSeq + dbEndPos + 1, dbFragLen);
                        if (fragment.size() + query.size() >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Sequence too long in query id: " << queryId << ". "
                                    "Max length allowed would is " << par.maxSeqLen << "\n";
                            break;
                        }
                        //update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                        queryCouldBeExtendedRight = true;
                        query += fragment;
                        rightQueryOffset += dbFragLen;

                    } else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                        if (queryCouldBeExtendedLeft == true) {
                            float alnLen = qEndPos - qStartPos;
                            float scorePerCol = static_cast<float>(score) / (alnLen+0.5);
                            besttHitToExtend.score = static_cast<int>(scorePerCol*100);
                            tmpAlignments.push_back(besttHitToExtend);
                            continue;
                        }
                        std::string fragment = std::string(targetSeq, dbStartPos); // +1 get not aligned element
                        if (fragment.size() + query.size() >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Sequence too long in query id: " << queryId << ". "
                                    "Max length allowed would is " << par.maxSeqLen << "\n";
                            break;
                        }
                        // update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                        queryCouldBeExtendedLeft = true;
                        query = fragment + query;
                        leftQueryOffset += dbStartPos;
                    }

                }
                if (queryCouldBeExtendedRight || queryCouldBeExtendedLeft){
                    queryCouldBeExtended = true;
                }
                alignments.clear();
                querySeqLen = query.size();
                querySeq = (char *) query.c_str();
                for(size_t alnIdx = 0; alnIdx < tmpAlignments.size(); alnIdx++){
                    int idCnt = 0;
                    int qStartPos = tmpAlignments[alnIdx].qStartPos;
                    int qEndPos = tmpAlignments[alnIdx].qEndPos;
                    int dbStartPos = tmpAlignments[alnIdx].dbStartPos;
                    int diagonal = (leftQueryOffset + besttHitToExtend.qStartPos) - besttHitToExtend.dbStartPos;
                    int dist = std::max(abs(diagonal), 0);
                    if (diagonal >= 0) {
                        qStartPos+=dist;
                        qEndPos+=dist;
                    }else{
                        dbStartPos+=dist;
                    }
                    unsigned int targetId = sequenceDbr->getId(tmpAlignments[alnIdx].dbKey);
                    char *targetSeq = sequenceDbr->getData(targetId);
                    for(int i = qStartPos; i <= qEndPos; i++){
                        int targetRes = static_cast<int>(targetSeq[dbStartPos+(i-qStartPos)]);
                        int queryRes = static_cast<int>(querySeq[i]);
                        idCnt += (queryRes == targetRes) ? 1 : 0;
                    }
                    float seqId =  static_cast<float>(idCnt) / (static_cast<float>(qEndPos) - static_cast<float>(qStartPos) + 0.5);
                    tmpAlignments[alnIdx].seqId = seqId;
                    if(seqId >= par.seqIdThr){
                        alignments.push_back(tmpAlignments[alnIdx]);
                    }
                }
            }
            if (queryCouldBeExtended == true) {
                query.push_back('\n');
                __sync_or_and_fetch(&wasExtended[id], static_cast<unsigned char>(0x20));
                resultWriter.writeData(query.c_str(), query.size(), queryId, thread_idx);
            }
        }
    } // end parallel

// add sequences that are not yet assembled
#pragma omp parallel for schedule(dynamic, 10000)
    for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        //   bool couldExtend =  (wasExtended[id] & 0x10);
        bool isNotContig =  !(wasExtended[id] & 0x20);
//        bool wasNotUsed =  !(wasExtended[id] & 0x40);
//        bool wasNotExtended =  !(wasExtended[id] & 0x80);
        //    bool wasUsed    =  (wasExtended[id] & 0x40);
        //if(isNotContig && wasNotExtended ){
        if (isNotContig){
            char *querySeqData = sequenceDbr->getData(id);
            unsigned int queryLen = sequenceDbr->getSeqLens(id) - 1; //skip null byte
            resultWriter.writeData(querySeqData, queryLen, sequenceDbr->getDbKey(id), thread_idx);
        }
    }

    // cleanup

    resultWriter.close(sequenceDbr->getDbtype());
    alnReader->close();
    delete [] wasExtended;
    delete alnReader;
    delete [] fastMatrix.matrix;
    delete [] fastMatrix.matrixData;
    sequenceDbr->close();
    delete sequenceDbr;
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}

int assembleresult(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, 3);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Compute assembly.\n";
    return doassembly(par);
}

