#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "tmalign/basic_fun.h"
#include "MultimerUtil.h"
#include "LDDT.h"
#include <map>
#include <set>

#ifdef OPENMP
#include <omp.h>
#endif

struct Complex {
    int complexId;
    unsigned int nChain;
    unsigned int complexLength;
    std::string complexName;
    std::vector<unsigned int> chainLengths;
    std::vector<unsigned int> chainKeys;

    // Coordinate16 Coords;

    Complex() : complexId(0), nChain(0), complexLength(0), complexName("") {}
    ~Complex() {
        chainKeys.clear();
    }
};

struct AlignedCoordinate {
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    AlignedCoordinate() {}
    AlignedCoordinate(size_t size) {
        x.resize(size);
        y.resize(size);
        z.resize(size);
    }
    ~AlignedCoordinate() {
        x.clear();
        y.clear();
        z.clear();
    }
};

unsigned int adjustAlnLen(unsigned int qcov, unsigned int tcov, int covMode) {
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return (qcov+tcov)/2;
        case Parameters::COV_MODE_TARGET:
            return qcov;
        case Parameters::COV_MODE_QUERY:
            return tcov;
        default:
            return 0;
    }
}

class ComplexFilterCriteria {
public:
    unsigned int targetComplexId;

    // per complex
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    float qCov;
    float tCov;
    float interfaceLddt;
    double qTm;
    double tTm;
    double avgTm;
    float t[3];
    float u[3][3];

    // per chain : criteria for chainTmThr & lddtThr
    std::vector<unsigned int> qAlnChainKeys;
    std::vector<unsigned int> tAlnChainKeys;
    std::vector<AlignedCoordinate> qAlnChains;
    std::vector<AlignedCoordinate> tAlnChains;

    std::vector<double> qAlnChainTms;
    std::vector<double> tAlnChainTms;

    ComplexFilterCriteria() {}
    ComplexFilterCriteria(
        unsigned int targetComplexId, double qTm, double tTm, float tstring[3], float ustring[3][3]
    ) :
        targetComplexId(targetComplexId), qTotalAlnLen(0), tTotalAlnLen(0), qCov(0), tCov(0), interfaceLddt(0), qTm(qTm), tTm(tTm), avgTm(0) {
        std::copy(tstring, tstring + 3, t);
        for (int i = 0; i < 3; i++) {
            std::copy(ustring[i], ustring[i] + 3, u[i]);
        }
    }
    ~ComplexFilterCriteria() {
        qAlnChainTms.clear();
        tAlnChainTms.clear();
        qAlnChainKeys.clear();
        tAlnChainKeys.clear();
        qAlnChains.clear();
        tAlnChains.clear();
    }

    bool hasTm(float TmThr, int covMode) {
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                return ((qTm>= TmThr) && (tTm >= TmThr));
            case Parameters::COV_MODE_TARGET:   
                return (tTm >= TmThr);
            case Parameters::COV_MODE_QUERY:
                return (qTm >= TmThr);
            default:
                return true;
        }
    }

    bool hasChainTm(float chainTmThr, int covMode, unsigned int qChainNum, unsigned int tChainNum) {
        if (qAlnChainTms.size()<std::min(qChainNum, tChainNum)) {
            return false;
        }
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                for (size_t i = 0; i < qAlnChainTms.size(); i++) {
                    if (qAlnChainTms[i] < chainTmThr || tAlnChainTms[i] < chainTmThr) {
                        return false;
                    }
                }
                break;
            case Parameters::COV_MODE_TARGET:
                for (size_t i = 0; i < tAlnChainTms.size(); i++) {
                    if (tAlnChainTms[i] < chainTmThr) {
                        return false;
                    }
                }
                break;
            case Parameters::COV_MODE_QUERY:
                for (size_t i = 0; i < qAlnChainTms.size(); i++) {
                    if (qAlnChainTms[i] < chainTmThr) {
                        return false;
                    }
                }
                break;
            default:
                return true;
        }
        return true;
    }

    bool hasChainNum(int covMode, unsigned int qChainNum, unsigned int tChainNum) {
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                if (qChainNum != tChainNum) {
                    return false;
                }
                break;
            default:
                return true;
        }
        return true;
    }

    // void calculateAvgTm(int covMode){
    //     switch (covMode) {
    //         case Parameters::COV_MODE_BIDIRECTIONAL:
    //             avgTm = ( qTm + tTm ) / 2 ;
    //             break;
    //         case Parameters::COV_MODE_TARGET:
    //             avgTm = tTm ;
    //             break;
    //         case Parameters::COV_MODE_QUERY:
    //             avgTm = qTm ;
    //             break;
    //         default :
    //             avgTm = ( qTm + tTm ) / 2 ;
    //     }
    // }

    bool hasInterfaceLDDT(float iLddtThr, unsigned int qChainNum, unsigned int tChainNum) {
        if (qAlnChainTms.size()<std::min(qChainNum, tChainNum)) {
            return false;
        }
        return(interfaceLddt >= iLddtThr);
    }
    bool satisfy(int covMode, float covThr, float TmThr, float chainTmThr, float iLddtThr, size_t qChainNum, size_t tChainNum ) {
        const bool covOK = covThr ? Util::hasCoverage(covThr, covMode, qCov, tCov) : true;
        const bool TmOK = TmThr ? hasTm(TmThr, covMode) : true;
        const bool chainTmOK = chainTmThr ? hasChainTm(chainTmThr, covMode, qChainNum, tChainNum) : true;
        const bool chainNumOK = hasChainNum(covMode, qChainNum, tChainNum);
        const bool lddtOK = iLddtThr ? hasInterfaceLDDT(iLddtThr, qChainNum, tChainNum) : true;
        // calculateAvgTm(covMode);
        return (covOK && TmOK && chainTmOK && lddtOK && chainNumOK);
    }

    void updateAln(unsigned int qAlnLen, unsigned int tAlnLen) {
        qTotalAlnLen += qAlnLen;
        tTotalAlnLen += tAlnLen;
    }

    void updateChainTmScore(double qChainTm, double tChainTm) {
        qAlnChainTms.push_back(qChainTm);
        tAlnChainTms.push_back(tChainTm);
    }

    void fillChainAlignment(unsigned int qChainKey, unsigned int tChainKey, unsigned int alnLen, 
                            float *qdata, float *tdata, const std::string &cigar, int qStartPos, int tStartPos, int qLen, int tLen) {
        AlignedCoordinate qChain;
        AlignedCoordinate tChain;
        int qi = qStartPos;
        int ti = tStartPos;
        int mi = 0;
        std::string backtrace = Matcher::uncompressAlignment(cigar);

        qChain.x.resize(alnLen);
        qChain.y.resize(alnLen);
        qChain.z.resize(alnLen);
        tChain.x.resize(alnLen);
        tChain.y.resize(alnLen);
        tChain.z.resize(alnLen);        

        for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
            if (backtrace[btPos] == 'M') {
                qChain.x[mi] = qdata[qi];
                qChain.y[mi] = qdata[qLen + qi];
                qChain.z[mi] = qdata[2*qLen + qi];
                tChain.x[mi] = tdata[ti];
                tChain.y[mi] = tdata[tLen + ti];
                tChain.z[mi] = tdata[2*tLen + ti];
                qi++;
                ti++;
                mi++;
            }
            else if (backtrace[btPos] == 'I') {
                qi++;
            }
            else {
                ti++;
            }
        }
        qAlnChainKeys.push_back(qChainKey);
        tAlnChainKeys.push_back(tChainKey);
        qAlnChains.push_back(qChain);
        tAlnChains.push_back(tChain);
    }
    // void update(unsigned int qChainKey, unsigned int tChainKey, double qChainTm, double tChainTm) {
    //     this->qAlnChainTms.push_back(qChainTm);
    //     this->tAlnChainTms.push_back(tChainTm);
        
    //     this->qAlnChainKeys.push_back(qChainKey);
    //     this->tAlnChainKeys.push_back(tChainKey);
    // }

    void calcCov(unsigned int qLen, unsigned int tLen) {
        qCov = static_cast<float>(qTotalAlnLen) / static_cast<float>(qLen);
        tCov = static_cast<float>(tTotalAlnLen) / static_cast<float>(tLen);
    }

    void computeInterfaceLddt(float threshold = 8) {
        if (qAlnChains.size() == 1) {
            interfaceLddt = 1;
        }
        float t2 = threshold * threshold;
        std::vector<std::set<unsigned int>> qInterfacePos(qAlnChains.size()); // chainIdx, resIdx
        unsigned int intLen = 0;
        // Find and save interface Coordinates
        for (size_t chainIdx1 = 0; chainIdx1 < qAlnChains.size(); chainIdx1++) {
            for (size_t chainIdx2 = chainIdx1+1; chainIdx2 < qAlnChains.size(); chainIdx2++) {
                AlignedCoordinate qChain1 = qAlnChains[chainIdx1];
                AlignedCoordinate qChain2 = qAlnChains[chainIdx2];
                for (size_t resIdx1 = 0; resIdx1 < qChain1.x.size(); resIdx1++) {
                    for (size_t resIdx2 = 0; resIdx2 < qChain2.x.size(); resIdx2++) {
                        float dist = BasicFunction::dist(qChain1.x[resIdx1], qChain1.y[resIdx1], qChain1.z[resIdx1],
                                                         qChain2.x[resIdx2], qChain2.y[resIdx2], qChain2.z[resIdx2]);
                        if (dist < t2) {
                            if (qInterfacePos[chainIdx1].find(resIdx1) == qInterfacePos[chainIdx1].end()) {
                                qInterfacePos[chainIdx1].insert(resIdx1);
                                intLen++;
                            }
                            if (qInterfacePos[chainIdx2].find(resIdx2) == qInterfacePos[chainIdx2].end()) {
                                qInterfacePos[chainIdx2].insert(resIdx2);
                                intLen++;
                            }
                        }
                    }
                }
            }
        }

        if (intLen == 0) {
            return;
        }
        AlignedCoordinate qInterface(intLen);
        AlignedCoordinate tInterface(intLen);
        size_t idx = 0;
        for (size_t chainIdx = 0; chainIdx < qInterfacePos.size(); chainIdx++) {
            if (qInterfacePos[chainIdx].size() >= 4) {
                for (size_t resIdx: qInterfacePos[chainIdx]) {
                    qInterface.x[idx] = qAlnChains[chainIdx].x[resIdx];
                    qInterface.y[idx] = qAlnChains[chainIdx].y[resIdx];
                    qInterface.z[idx] = qAlnChains[chainIdx].z[resIdx];
                    tInterface.x[idx] = tAlnChains[chainIdx].x[resIdx];
                    tInterface.y[idx] = tAlnChains[chainIdx].y[resIdx];
                    tInterface.z[idx] = tAlnChains[chainIdx].z[resIdx];
                    idx++;
                }
            }
        }
        std::string bt(intLen, 'M');
        LDDTCalculator *lddtcalculator = NULL;
        lddtcalculator = new LDDTCalculator(intLen+1, intLen+1);
        lddtcalculator->initQuery(intLen, &qInterface.x[0], &qInterface.y[0], &qInterface.z[0]);
        LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(intLen, 0, 0, bt, &tInterface.x[0], &tInterface.y[0], &tInterface.z[0]);
        interfaceLddt = lddtres.avgLddtScore;
        delete lddtcalculator;
    }
};


char* filterToBuffer(ComplexFilterCriteria cmplfiltcrit, char* tmpBuff){
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.qCov, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.tCov, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.qTm, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.tTm, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.interfaceLddt, tmpBuff);    
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[0][0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[0][1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[0][2], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[1][0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[1][1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[1][2], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[2][0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[2][1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[2][2], tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.t[0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.t[1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.t[2], tmpBuff);
    *(tmpBuff-1) = '\n';
    return tmpBuff;
}

void fillUArr(const std::string &uString, float (&u)[3][3]) {
    std::string tmp;
    int i = 0;
    int j=0;
    const int ulen = static_cast<int>(uString.size());
    for (int k=0; k < ulen; k++) {
        if (k==ulen-1) {
            u[i][j] = std::stof(tmp);
        } else if (uString[k] == ',') {
            u[i][j] = std::stof(tmp);
            tmp.clear();
            j++;
        } else {
            tmp.push_back(uString[k]);
        }
        if (j == 3) {
            i++;
            j = 0;
        }
    }
}

void fillTArr(const std::string &tString, float (&t)[3]) {
    std::string tmp;
    int i = 0;
    const int tlen = static_cast<int>(tString.size());
    for (int k=0; k<tlen; k++) {
        if (k ==tlen-1) {
            t[i] = std::stof(tmp);
        } else if (tString[k] == ',') {
            t[i] = std::stof(tmp);
            tmp.clear();
            i++;
        } else {
            tmp.push_back(tString[k]);
        }
    }
}

unsigned int cigarToAlignedLength(const std::string &cigar) {
    std::string backtrace = Matcher::uncompressAlignment(cigar);
    unsigned int alni = 0;
    for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
        if (backtrace[btPos] == 'M') {
            alni++;
        }
    }
    return alni;
}

double computeChainTmScore(AlignedCoordinate &qchain, AlignedCoordinate &tchain, float t[3], float u[3][3], int tLen) {
    unsigned int alnLen = qchain.x.size();
    double tmscore = 0;
    float d0 = 1.24*(cbrt(tLen-15)) -1.8;
    float d02 = d0*d0;

    Coordinates tmt(alnLen);
    // BasicFunction::do_rotation(tchain.x, tchain.y, tchain.z, tmt, alnLen, t, u);
    for (unsigned int k=0; k<alnLen; k++) {
        float tmx, tmy, tmz;
        BasicFunction::transform(t, u, tchain.x[k], tchain.y[k], tchain.z[k], tmx, tmy, tmz);
        // double di = BasicFunction::dist(qchain.x[k], qchain.y[k], qchain.z[k], tmt.x[k], tmt.y[k], tmt.z[k]);
        double di = BasicFunction::dist(qchain.x[k], qchain.y[k], qchain.z[k], tmx, tmy, tmz);
        tmscore += 1/(1+di/d02);
    }
    return tmscore;
}

void getComplexResidueLength( IndexReader *Dbr, std::vector<Complex> &complexes) {
    for (size_t complexIdx = 0; complexIdx < complexes.size(); complexIdx++) {
        Complex *complex = &complexes[complexIdx];
        std::vector<unsigned int> &chainKeys = complex->chainKeys;
        if (chainKeys.empty()) {
            continue;
        }
        unsigned int cmpllen = 0;
        for (auto chainKey: chainKeys) {
            size_t id = Dbr->sequenceReader->getId(chainKey);
            if (id == NOT_AVAILABLE_CHAIN_KEY) {
                break;
            }
            unsigned int reslen = Dbr->sequenceReader->getSeqLen(id);
            complex->chainLengths.push_back(reslen);
            cmpllen += Dbr->sequenceReader->getSeqLen(id);
        }
        complex->complexLength = cmpllen;
    }
}

static void getlookupInfo(
        IndexReader* dbr,
        const std::string &file,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::vector<Complex> &complexes,
        std::map<unsigned int, unsigned int> &complexIdtoIdx
) {
    if (file.length() == 0) {
        return;
    }
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];

    int prevComplexId =  -1;
    int nComplex = 0;
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        auto chainKey = Util::fast_atoi<int>(entry[0]);
        unsigned int chainDbId = dbr->sequenceReader->getId(chainKey);
        if (chainDbId != NOT_AVAILABLE_CHAIN_KEY) {
            auto complexId = Util::fast_atoi<int>(entry[2]);
            chainKeyToComplexIdLookup.emplace(chainKey, complexId);
            std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
            size_t lastUnderscoreIndex = chainName.find_last_of('_');
            std::string complexName = chainName.substr(0, lastUnderscoreIndex);

            if (complexId != prevComplexId) {
                
                Complex complex;
                complex.complexId = complexId;
                complex.complexName = complexName;
                complexIdtoIdx.emplace(complexId, nComplex);
                complexes.emplace_back(complex);

                prevComplexId = complexId;
                nComplex++;
            }
            complexes.back().chainKeys.emplace_back(chainKey);
            complexes.back().nChain++;
        }
        data = Util::skipLine(data);
    }
    lookupDB.close();
}

int filtermultimer(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);

    IndexReader* qDbr = NULL;
    qDbr = new IndexReader(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qStructDbr.open(DBReader<unsigned int>::NOSORT);

    IndexReader* tDbr = NULL;
    DBReader<unsigned int> *tStructDbr = NULL;
    if (sameDB) {
        tDbr = qDbr;
        tStructDbr = &qStructDbr;
    }
    else{
        tDbr = new IndexReader(par.db2, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tStructDbr = new DBReader<unsigned int>((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(),
                                           par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tStructDbr->open(DBReader<unsigned int>::NOSORT);
    }
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX| DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t localThreads = 1;
#ifdef OPENMP
localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif
    const bool shouldCompress = (par.compressed == true);
    const int db4Type = Parameters::DBTYPE_CLUSTER_RES;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, shouldCompress, db4Type);
    resultWriter.open();
    const int db5Type = Parameters::DBTYPE_GENERIC_DB;
    DBWriter resultWrite5((par.db4 + "_info").c_str(), (par.db4 + "_info.index").c_str(), par.threads, shouldCompress, db5Type);
    resultWrite5.open();

    std::string qLookupFile = par.db1 + ".lookup";
    std::string tLookupFile = par.db2 + ".lookup";
    
    std::vector<Complex> qComplexes, tComplexes;
    std::map<unsigned int, unsigned int> qComplexIdToIdx, tComplexIdToIdx;
    chainKeyToComplexId_t qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;

    getlookupInfo(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexes, qComplexIdToIdx);
    getComplexResidueLength(qDbr, qComplexes);
    Debug::Progress progress(qComplexes.size());

    if (sameDB) {
        tChainKeyToComplexIdMap = qChainKeyToComplexIdMap;
        tComplexes = qComplexes;
        tComplexIdToIdx = qComplexIdToIdx;
    } else {
        getlookupInfo(tDbr, tLookupFile, tChainKeyToComplexIdMap, tComplexes, tComplexIdToIdx);
        getComplexResidueLength(tDbr, tComplexes);
    }
    // std::vector<unsigned int> qComplexOrder(qComplexes.size());
    // for (size_t qComplexIdx = 0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
    //         qComplexOrder[qComplexIdx] = qComplexIdx;
    // }
    // std::sort(qComplexOrder.begin(), qComplexOrder.end(), [&qComplexes](unsigned int lhs, unsigned int rhs) {
    //           return qComplexes[lhs].chainKeys.size() > qComplexes[rhs].chainKeys.size();
    //       });
    
#pragma omp parallel num_threads(localThreads) 
    {   
        char buffer[32];
        char buffer2[4096];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t result;
        std::map<unsigned int, ComplexFilterCriteria> localComplexMap;
        // std::vector< ComplexFilterCriteria> localComplexVector;
        std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId;
        // std::vector<unsigned int> cmpltargetIds;
        // std::vector<double> targetIdBestTm;
        std::vector<unsigned int> selectedAssIDs;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        
        Matcher::result_t res;   
#pragma omp for schedule(dynamic, 1) 
        for (size_t qComplexIdx = 0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
        // for (size_t qComplexIdx : qComplexOrder) {
            progress.updateProgress();
            Complex qComplex = qComplexes[qComplexIdx];
            unsigned int qComplexId = qComplex.complexId;
            std::vector<unsigned int> qChainKeys = qComplex.chainKeys;
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainAlnId = alnDbr.getId(qChainKey);
                unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                // Handling monomer as singleton
                if (qChainAlnId == NOT_AVAILABLE_CHAIN_KEY) {
                    break;
                }
                
                char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
                size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
                size_t qChainLen = qDbr->sequenceReader->getSeqLen(qChainDbId);
                float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);
                
                char *data = alnDbr.getData(qChainAlnId, thread_idx);
                while (*data != '\0' ) {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                
                    if (!retComplex.isValid) {
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }

                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainKey = res.dbKey;
                    unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                    unsigned int tComplexId = tChainKeyToComplexIdMap.at(tChainKey);
                    unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                    std::vector<unsigned int> tChainKeys = tComplexes[tComplexIdx].chainKeys;
                    //if target is monomer, but user doesn't want, continue
                    unsigned int tChainAlnId = alnDbr.getId(tChainKey);
                    if (tChainAlnId == NOT_AVAILABLE_CHAIN_KEY) {
                        continue;
                    }
                    float u[3][3];
                    float t[3];
                    fillUArr(retComplex.uString, u);
                    fillTArr(retComplex.tString, t);
                    unsigned int qalnlen = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    unsigned int talnlen = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                    // if (localComplexVector.size() <= assId) {
                    //     ComplexFilterCriteria cmplfiltcrit(tComplexId, retComplex.qTmScore, retComplex.tTmScore, t, u);
                    //     size_t subt = assId - localComplexVector.size();
                    //     for (size_t sub=0; sub < subt; sub ++) {
                    //         localComplexVector.push_back(cmplfiltcrit);
                    //     }
                    //     localComplexVector.push_back(cmplfiltcrit);
                    // }
                    if (localComplexMap.find(assId) == localComplexMap.end()) {
                        ComplexFilterCriteria cmplfiltcrit(tComplexId, retComplex.qTmScore, retComplex.tTmScore, t, u);
                        localComplexMap[assId] = cmplfiltcrit;
                    }
                    ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);
                    // ComplexFilterCriteria &cmplfiltcrit = localComplexVector.at(assId);
                    cmplfiltcrit.updateAln(qalnlen, talnlen);
    
                    // save Aligned coordinatese if needed : chainTmThr & lddtThr
                    if (par.filtChainTmThr > 0.0f || par.filtInterfaceLddtThr > 0.0f) {
                        char *tcadata = tStructDbr->getData(tChainDbId, thread_idx);
                        size_t tCaLength = tStructDbr->getEntryLen(tChainDbId);
                        float* tdata = tcoords.read(tcadata, res.dbLen, tCaLength);

                        unsigned int alnLen = cigarToAlignedLength(res.backtrace);
                        cmplfiltcrit.fillChainAlignment(qChainKey, tChainKey, alnLen, qdata, tdata, res.backtrace, res.qStartPos, res.dbStartPos, res.qLen, res.dbLen);
                        // if (par.filtChainTmThr > 0.0f) {
                        double chainTm = computeChainTmScore(cmplfiltcrit.qAlnChains.back(), cmplfiltcrit.tAlnChains.back(), t, u, res.dbLen);
                        cmplfiltcrit.updateChainTmScore(chainTm / res.qLen, chainTm / res.dbLen);
                        // }
                    }
                } // while end
            }
            
            // Filter the target complexes and get the best alignment
            // for (unsigned int assId = 0; assId < localComplexVector.size(); assId++) {
            for (auto& assId_res : localComplexMap) {
                unsigned int tComplexId  = assId_res.second.targetComplexId;
                // unsigned int tComplexId  = localComplexVector.at(assId).targetComplexId;
                
                unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                Complex  tComplex = tComplexes[tComplexIdx];

                ComplexFilterCriteria &cmplfiltcrit = assId_res.second;
                // ComplexFilterCriteria &cmplfiltcrit = localComplexVector.at(assId);
                cmplfiltcrit.calcCov(qComplex.complexLength, tComplex.complexLength);

                if (par.filtInterfaceLddtThr > 0.0) {
                    cmplfiltcrit.computeInterfaceLddt();
                }

                // Check if the criteria are met
                if (!(cmplfiltcrit.satisfy(par.covMode, par.covThr, par.filtMultimerTmThr, par.filtChainTmThr, par.filtInterfaceLddtThr, qComplex.nChain, tComplex.nChain))) {
                    continue;
                }
                unsigned int alnlen = adjustAlnLen(cmplfiltcrit.qTotalAlnLen, cmplfiltcrit.tTotalAlnLen, par.covMode);
                // Get the best alignement per each target complex   
                if (cmplIdToBestAssId.find(tComplexId) == cmplIdToBestAssId.end()) {
                    cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    // cmplIdToBestAssId[tComplexId] = {static_cast<double>(assId_res.first), cmplfiltcrit.avgTm};
                    // cmplIdToBestAssId[tComplexId] = {static_cast<double>(assId), cmplfiltcrit.avgTm};
                } else {
                    if (alnlen > cmplIdToBestAssId.at(tComplexId)[1]) {
                        cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    }
                    // if (cmplfiltcrit.avgTm > cmplIdToBestAssId.at(tComplexId)[1]) {
                    //     cmplIdToBestAssId[tComplexId] = {static_cast<double>(assId_res.first), cmplfiltcrit.avgTm};
                    //     // cmplIdToBestAssId[tComplexId] = {static_cast<double>(assId_res), cmplfiltcrit.avgTm};
                    // }
                }
                
                // unsigned int targetindex;
                // auto it = std::find(cmpltargetIds.begin(), cmpltargetIds.end(), tComplexId);
                // if ( it == cmpltargetIds.end()) {
                //     cmpltargetIds.push_back(tComplexId);
                //     selectedAssIDs.push_back(assId);
                //     targetIdBestTm.push_back(cmplfiltcrit.avgTm);
                // } else {
                //     targetindex = std::distance(cmpltargetIds.begin(), it);
                //     if (cmplfiltcrit.avgTm > targetIdBestTm[targetindex]) {
                //         targetIdBestTm[targetindex] = cmplfiltcrit.avgTm;
                //         selectedAssIDs[targetindex] = assId;
                //     }
                // }

            }

            for (const auto& pair : cmplIdToBestAssId) {
                selectedAssIDs.push_back(pair.second[0]);
            }
            resultWrite5.writeStart(thread_idx);
            for (unsigned int assIdidx = 0; assIdidx < selectedAssIDs.size(); assIdidx++) {
                unsigned int assId = selectedAssIDs.at(assIdidx);
                ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);
                // ComplexFilterCriteria &cmplfiltcrit = localComplexVector.at(assId);
                unsigned int tComplexId = cmplfiltcrit.targetComplexId;
                unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                Complex tComplex = tComplexes.at(tComplexIdx);
                
                char *outpos = Itoa::u32toa_sse2(tComplexId, buffer);
                result.append(buffer, (outpos - buffer - 1));
                result.push_back('\n');

                char * tmpBuff = Itoa::u32toa_sse2(tComplexId, buffer2);
                tmpBuff = filterToBuffer(cmplfiltcrit, tmpBuff);
                resultWrite5.writeAdd(buffer2, tmpBuff - buffer2, thread_idx);
            }
            if (selectedAssIDs.size() == 0) {
                float t[3];
                float u[3][3];
                for (int i=0; i < 3; i++) {
                    t[i] = 0.0;
                }
                for (int i=0; i < 3; i++) {
                    for (int j=0; j < 3; j++) {
                        u[i][j] = 0.0;
                    }
                }
                ComplexFilterCriteria cmplfiltcrit(qComplexId, 1.0, 1.0, t, u);
                cmplfiltcrit.qCov = 1.0;
                cmplfiltcrit.tCov = 1.0;
                cmplfiltcrit.interfaceLddt = 1.0;
                resultWrite5.writeStart(thread_idx);
                char * tmpBuff = Itoa::u32toa_sse2(qComplexId, buffer2);
                tmpBuff = filterToBuffer(cmplfiltcrit, tmpBuff);
                resultWrite5.writeAdd(buffer2, tmpBuff - buffer2, thread_idx);

                char *outpos = Itoa::u32toa_sse2(qComplexId, buffer);
                result.append(buffer, (outpos - buffer - 1));
                result.push_back('\n');
            }
            // if (qComplexId == 1) {
            //     Debug(Debug::WARNING)<< "hi\n";
            // }
            resultWriter.writeData(result.c_str(), result.length(), qComplexId, thread_idx);
            resultWrite5.writeEnd(qComplexId, thread_idx);
            result.clear();
            localComplexMap.clear();
            cmplIdToBestAssId.clear();
            selectedAssIDs.clear();        
            // localComplexVector.clear();
            // cmpltargetIds.clear();
            // targetIdBestTm.clear();
        } // for end
    } // MP end
    
    resultWriter.close(true);
    resultWrite5.close(true);
    qStructDbr.close();
    alnDbr.close();
    delete qDbr;
    if (sameDB == false) {
        delete tDbr;
        delete tStructDbr;
    }
    qChainKeyToComplexIdMap.clear();
    tChainKeyToComplexIdMap.clear();
    qComplexes.clear();
    tComplexes.clear();
    return EXIT_SUCCESS;
}