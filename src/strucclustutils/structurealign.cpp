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
#include "TMaligner.h"
#include "Coordinate16.h"
#include "LDDT.h"

#ifdef OPENMP
#include <omp.h>
#endif

// need for sorting the results
static bool compareHitsByStructureBits(const Matcher::result_t &first, const Matcher::result_t &second) {
    if (first.score != second.score) {
        return first.score > second.score;
    }
    if (first.dbLen != second.dbLen) {
        return first.dbLen < second.dbLen;
    }
    return first.dbKey < second.dbKey;
}


static void structureAlignDefault(LocalParameters & par) {
    par.compBiasCorrectionScale = 0.5;
    par.alignmentType = LocalParameters::ALIGNMENT_TYPE_3DI_AA;
}

int alignStructure(StructureSmithWaterman & structureSmithWaterman,
                   StructureSmithWaterman & reverseStructureSmithWaterman,
                   Sequence & tSeqAA, Sequence & tSeq3Di,
                   unsigned int querySeqLen, unsigned int targetSeqLen,
                   EvalueNeuralNet & evaluer, std::pair<double, double> muLambda,
                   Matcher::result_t & res, std::string & backtrace,
                   Parameters & par) {

    float seqId = 0.0;
    backtrace.clear();
    // align only score and end pos
    StructureSmithWaterman::s_align align = structureSmithWaterman.alignScoreEndPos<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence, tSeq3Di.numSequence, targetSeqLen, par.gapOpen.values.aminoacid(),
                                                                                    par.gapExtend.values.aminoacid(), querySeqLen / 2);
    bool hasLowerCoverage = !(Util::hasCoverage(par.covThr, par.covMode, align.qCov, align.tCov));
    if(hasLowerCoverage){
        return -1;
    }
    // we can already stop if this e-value isn't good enough, it wont be any better in the next step
    align.evalue = evaluer.computeEvalueCorr(align.score1, muLambda.first, muLambda.second);
    bool hasLowerEvalue = align.evalue > par.evalThr;
    if(hasLowerEvalue){
        return -1;
    }

    StructureSmithWaterman::s_align revAlign;
    if(structureSmithWaterman.isProfileSearch()){
        revAlign.score1 = 0;
    } else {
        revAlign = reverseStructureSmithWaterman.alignScoreEndPos<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence, tSeq3Di.numSequence,
                                                                  targetSeqLen, par.gapOpen.values.aminoacid(),
                                                                  par.gapExtend.values.aminoacid(), querySeqLen / 2);
    }
    int32_t score = static_cast<int32_t>(align.score1) - static_cast<int32_t>(revAlign.score1);
    align.evalue = evaluer.computeEvalueCorr(score, muLambda.first, muLambda.second);
    hasLowerEvalue = align.evalue > par.evalThr;
    if (hasLowerEvalue) {
        return -1;
    }
    if (structureSmithWaterman.isProfileSearch()) {
        align = structureSmithWaterman.alignStartPosBacktrace<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence,
                                                                                               tSeq3Di.numSequence,
                                                                                               targetSeqLen,
                                                                                               par.gapOpen.values.aminoacid(),
                                                                                               par.gapExtend.values.aminoacid(),
                                                                                               par.alignmentMode,
                                                                                               backtrace, align,
                                                                                               par.covMode, par.covThr,
                                                                                               querySeqLen / 2);
    }else{
        align = structureSmithWaterman.alignStartPosBacktraceBlock(tSeqAA.numSequence, tSeq3Di.numSequence, targetSeqLen, par.gapOpen.values.aminoacid(),
                                                                   par.gapExtend.values.aminoacid(), backtrace, align);
    }

    unsigned int alnLength = Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1, align.dbStartPos1, align.dbEndPos1);
    if(backtrace.size() > 0){
        alnLength = backtrace.size();
        seqId = Util::computeSeqId(par.seqIdMode, align.identicalAACnt, querySeqLen, targetSeqLen, alnLength);
    }
    align.score1 = score;
    res = Matcher::result_t(tSeqAA.getDbKey(), align.score1, align.qCov, align.tCov, seqId, align.evalue, alnLength,
                            align.qStartPos1, align.qEndPos1, querySeqLen, align.dbStartPos1, align.dbEndPos1, targetSeqLen, backtrace);
    return 0;
}


int computeAlternativeAlignment(StructureSmithWaterman & structureSmithWaterman,
                                StructureSmithWaterman & reverseStructureSmithWaterman,
                                Sequence & tSeqAA, Sequence & tSeq3Di,
                                unsigned int querySeqLen, unsigned int targetSeqLen,
                                EvalueNeuralNet & evaluer, std::pair<double, double> muLambda,
                                Matcher::result_t & result, Matcher::result_t & altRes,
                                std::string & backtrace, Parameters & par) {
    const unsigned char xAAIndex = tSeqAA.subMat->aa2num[static_cast<int>('X')];
    const unsigned char x3DiIndex = tSeq3Di.subMat->aa2num[static_cast<int>('X')];
    for (int pos = result.dbStartPos; pos < result.dbEndPos; ++pos) {
        tSeqAA.numSequence[pos] = xAAIndex;
        tSeq3Di.numSequence[pos] = x3DiIndex;
    }
    if (alignStructure(structureSmithWaterman, reverseStructureSmithWaterman,
                       tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                       evaluer, muLambda, altRes, backtrace, par) == -1) {
        return -1;
    }
    if (Alignment::checkCriteria(altRes, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
        return 0;
    } else {
        return -1;
    }
}


int structurealign(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    structureAlignDefault(par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbrAA(par.db1, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    IndexReader qdbr3Di(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *t3DiDbr = NULL;
    IndexReader *tAADbr = NULL;
    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;

    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &qdbr3Di;
        tAADbr = &qdbrAA;
    } else {
        tAADbr = new IndexReader(par.db2, par.threads,
                                 alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                                 (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

        std::string t3DiDbrName =  StructureUtil::getIndexWithSuffix(par.db2, "_ss");
        bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()),
                                                  Parameters::DBTYPE_INDEX_DB);

        t3DiDbr = new IndexReader(is3DiIdx ? t3DiDbrName : par.db2, par.threads,
                                  alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                                  (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                                  alignmentIsExtended ? "_seq_ss" : "_ss");
    }

    bool db1CaExist = FileUtil::fileExists((par.db1 + "_ca.dbtype").c_str());
    bool db2CaExist = FileUtil::fileExists((par.db2 + "_ca.dbtype").c_str());
    if(Parameters::isEqualDbtype(tAADbr->getDbtype(), Parameters::DBTYPE_INDEX_DB)){
        db2CaExist = true;
    }
    if(par.sortByStructureBits) {
        bool disableStructureBits = false;
        if(db1CaExist == false || db2CaExist == false){
            Debug(Debug::WARNING) << "Cannot find " << FileUtil::baseName(par.db1) << " C-alpha or " << FileUtil::baseName(par.db2) << " C-alpha database\n";
            disableStructureBits = true;
        }
        if(par.alignmentMode == 1 || par.alignmentMode == 2){
            Debug(Debug::WARNING) << "Cannot use --sort-by-structure-bits 1 with --alignment-mode 1 or 2\n";
            disableStructureBits = true;
        }
        if(disableStructureBits){
            Debug(Debug::WARNING) << "Disabling --sort-by-structure-bits\n";
            Debug(Debug::WARNING) << "This impacts the final score and ranking of hits, but not E-values themselves. Ranking alterations primarily occur for E-values < 10^-1.\n";
            par.sortByStructureBits = false;
        }
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype =  Parameters::DBTYPE_ALIGNMENT_RES;
    if(alignmentIsExtended){
        dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  dbtype);
    dbw.open();

    bool needTMaligner = (par.tmScoreThr > 0);
    bool needLDDT = (par.lddtThr > 0);
    if(par.sortByStructureBits){
        needLDDT = true;
        needTMaligner = true;
    }
    bool needCalpha = (needTMaligner || needLDDT);
    IndexReader *qcadbr = NULL;
    IndexReader *tcadbr = NULL;
    if(needCalpha){
        qcadbr = new IndexReader(
                par.db1,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca");
        if (sameDB) {
            tcadbr = qcadbr;
        } else {
             tcadbr = new IndexReader(
                    par.db2,
                    par.threads,
                    alignmentIsExtended ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2) :
                                           IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                    touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                    DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                    alignmentIsExtended ? "_seq_ca" : "_ca"
            );
        }
    }

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
    float aaFactor = (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA) ? 1.4 : 0.0;
    SubstitutionMatrix subMatAA(blosum.c_str(), aaFactor, par.scoreBias);
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

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        EvalueNeuralNet evaluer(tAADbr->sequenceReader->getAminoAcidDBSize(), &subMat3Di);
        std::vector<Matcher::result_t> alignmentResult;
        StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMatAA, &subMat3Di);
        StructureSmithWaterman reverseStructureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMatAA, &subMat3Di);
        TMaligner *tmaligner = NULL;
        if(needTMaligner) {
            tmaligner = new TMaligner(
                    std::max(qdbr3Di.sequenceReader->getMaxSeqLen() + 1, t3DiDbr->sequenceReader->getMaxSeqLen() + 1), false, true);
        }
        LDDTCalculator *lddtcalculator = NULL;
        if(needLDDT) {
            lddtcalculator = new LDDTCalculator(qdbr3Di.sequenceReader->getMaxSeqLen() + 1,  t3DiDbr->sequenceReader->getMaxSeqLen() + 1);
        }
        Sequence qSeqAA(par.maxSeqLen, qdbrAA.getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, qdbr3Di.getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        std::string backtrace;
        char buffer[1024+32768];
        std::string resultBuffer;

        Coordinate16 qcoords;
        Coordinate16 tcoords;

        TMaligner::TMscoreResult tmres;
        LDDTCalculator::LDDTScoreResult lddtres;
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
                qSeq3Di.mapSequence(id, queryKey, querySeq3Di, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);
                if(needCalpha){
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
                std::pair<double, double> muLambda = evaluer.predictMuLambda(qSeq3Di.numSequence, qSeq3Di.L);
                structureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
                qSeq3Di.reverse();
                qSeqAA.reverse();
                reverseStructureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = t3DiDbr->sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                    char * targetSeq3Di = t3DiDbr->sequenceReader->getData(targetId, thread_idx);
                    char * targetSeqAA = tAADbr->sequenceReader->getData(targetId, thread_idx);
                    const int targetSeqLen = static_cast<int>(t3DiDbr->sequenceReader->getSeqLen(targetId));

                    tSeq3Di.mapSequence(targetId, dbKey, targetSeq3Di, targetSeqLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);
                    if(Util::canBeCovered(par.covThr, par.covMode, qSeq3Di.L, targetSeqLen) == false){
                        rejected++;
                        continue;
                    }
                    Matcher::result_t res;
                    if(alignStructure(structureSmithWaterman, reverseStructureSmithWaterman,
                                      tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                                      evaluer, muLambda, res, backtrace, par) == -1){
                        rejected++;
                        continue;
                    }

                    if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        if(needCalpha) {
                            size_t tId = tcadbr->sequenceReader->getId(res.dbKey);
                            char *tcadata = tcadbr->sequenceReader->getData(tId, thread_idx);
                            size_t tCaLength = tcadbr->sequenceReader->getEntryLen(tId);
                            float* targetCaData = tcoords.read(tcadata, res.dbLen, tCaLength);
                            if(needTMaligner) {
                                tmres = tmaligner->computeTMscore(targetCaData,
                                                                  &targetCaData[res.dbLen],
                                                                  &targetCaData[res.dbLen +
                                                                                res.dbLen],
                                                                  res.dbLen,
                                                                  res.qStartPos,
                                                                  res.dbStartPos,
                                                                  res.backtrace,
                                                                  std::min(static_cast<unsigned int>(res.backtrace.size()), std::min(res.dbLen, res.qLen)));
                                if (tmres.tmscore < par.tmScoreThr) {
                                    continue;
                                }
                            }
                            if(needLDDT){
                                lddtres = lddtcalculator->computeLDDTScore(res.dbLen, res.qStartPos, res.dbStartPos,
                                                                           res.backtrace,
                                                                           targetCaData, &targetCaData[res.dbLen],
                                                                           &targetCaData[res.dbLen+res.dbLen]);

                                if(lddtres.avgLddtScore < par.lddtThr){
                                    continue;
                                }
                                res.dbcov = lddtres.avgLddtScore;
                            }
                            if(par.sortByStructureBits && needTMaligner && needLDDT){
                                res.score = res.score * sqrt(lddtres.avgLddtScore * tmres.tmscore);
                            }
                        }


                        alignmentResult.emplace_back(res);
                        int altAli = par.altAlignment;
                        bool moreAltAli = true;
                        while(altAli && moreAltAli){
                            Matcher::result_t altRes;
                            if(computeAlternativeAlignment(structureSmithWaterman, reverseStructureSmithWaterman,
                                                           tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                                                           evaluer, muLambda, res, altRes,
                                                           backtrace, par) == -1) {
                                moreAltAli = false;
                                continue;
                            }
                            alignmentResult.push_back(altRes);
                            res = altRes;
                            altAli--;
                        }
                        passedNum++;
                        rejected = 0;
                    } else {
                        rejected++;
                    }
                }
            }


            if (alignmentResult.size() > 1) {
                if(par.sortByStructureBits) {
                    SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), compareHitsByStructureBits);
                } else {
                    SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), Matcher::compareHits);
                }
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

    dbw.close();
    resultReader.close();

    if(needCalpha){
        if (sameDB == false) {
            delete tcadbr;
        }
        delete qcadbr;
    }

    if (sameDB == false) {
        delete t3DiDbr;
        delete tAADbr;
    }

    return EXIT_SUCCESS;
}
