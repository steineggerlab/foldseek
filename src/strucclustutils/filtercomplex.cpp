#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "MemoryMapped.h"
#include "Coordinate16.h"
#include "tmalign/basic_fun.h"
#include "createcomplexreport.h"
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>


#ifdef OPENMP
#include <omp.h>
#endif


unsigned int adjustAlnLen(unsigned int qcov, unsigned int tcov, int covMode) {
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return (qcov+tcov)/2;
        case Parameters::COV_MODE_TARGET:
            return qcov;
        case Parameters::COV_MODE_QUERY:
            return tcov;
        case Parameters::COV_MODE_LENGTH_QUERY :
        case Parameters::COV_MODE_LENGTH_TARGET :
        case Parameters::COV_MODE_LENGTH_SHORTER :
            return 0;
    }
}

// bool hasChainTm(float chainTMThr, int covMode, std::vector<double> &qChainTmScores, std::vector<double> &tChainTmScores, unsigned int qChainNum, unsigned int tChainNum) {
//     if (chainTMThr > 0 ){
//         switch (covMode) {
//             case Parameters::COV_MODE_BIDIRECTIONAL:
//                 if (qChainTmScores.size()<std::min(qChainNum, tChainNum)){
//                     return false;
//                 }
//                 for (size_t i = 0; i < qChainTmScores.size(); i++) {
//                     if (qChainTmScores[i] < chainTMThr || tChainTmScores[i] < chainTMThr) {
//                         return false;
//                     }
//                 }
//                 break;
//             case Parameters::COV_MODE_TARGET:
//                 if (qChainTmScores.size()<std::min(qChainNum, tChainNum)){
//                     return false;
//                 }
//                 for (size_t i = 0; i < qChainTmScores.size(); i++) {
//                     if (tChainTmScores[i] < chainTMThr) {
//                         return false;
//                     }
//                 }
//                 break;
//             case Parameters::COV_MODE_QUERY:
//                 if (qChainTmScores.size()<std::min(qChainNum, tChainNum)){
//                     return false;
//                 }
//                 for (size_t i = 0; i < qChainTmScores.size(); i++) {
//                     if (qChainTmScores[i] < chainTMThr) {
//                         return false;
//                     }
//                 }
//                 break;
//             case Parameters::COV_MODE_LENGTH_QUERY :
//             case Parameters::COV_MODE_LENGTH_TARGET :
//             case Parameters::COV_MODE_LENGTH_SHORTER :
//                 break;
//         }
//     }
//     return true;
// }

int determineOctant(const float coord[3]) {
    if (coord[0] >= 0 && coord[1] >= 0 && coord[2] >= 0) return 0;
    if (coord[0] < 0 && coord[1] >= 0 && coord[2] >= 0) return 1;
    if (coord[0] < 0 && coord[1] < 0 && coord[2] >= 0) return 2;
    if (coord[0] >= 0 && coord[1] < 0 && coord[2] >= 0) return 3;
    if (coord[0] >= 0 && coord[1] >= 0 && coord[2] < 0) return 4;
    if (coord[0] < 0 && coord[1] >= 0 && coord[2] < 0) return 5;
    if (coord[0] < 0 && coord[1] < 0 && coord[2] < 0) return 6;
    if (coord[0] >= 0 && coord[1] < 0 && coord[2] < 0) return 7;
}

static bool isSameCoord(float qavgCoords[3], float tavgCoords[3]){
    return(determineOctant(qavgCoords) == determineOctant(tavgCoords)) ;
}

class ComplexFilterCriteria {
public:
    unsigned int dbKey;
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    float qCov;
    float tCov;
    double qTM;
    double tTM;
    bool sameCoord;
    std::vector<int> qOctaCount;
    std::vector<int> tOctaCount;
    std::vector<double> alignedQChainTmScores;
    std::vector<double> alignedTChainTmScores;
    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> tChainKeys;
    float refCoord[3];
    float t[3];
    float u[3][3];

    
    ComplexFilterCriteria() {}
    ComplexFilterCriteria(unsigned int dbKey, std::vector<unsigned int> &qChainKeys, std::vector<unsigned int> &tChainKeys, double qTM, double tTM, float tstring[3], float ustring[3][3]) :
                            dbKey(dbKey), qTM(qTM), tTM(tTM), qChainKeys(qChainKeys), tChainKeys(tChainKeys), qTotalAlnLen(0), tTotalAlnLen(0) {
                                std::copy(tstring, tstring + 3, t);
                                for (int i = 0; i < 3; i++) {
                                    std::copy(ustring[i], ustring[i] + 3, u[i]);
                                }
                                for (int i = 0; i < 8; i++) {
                                    qOctaCount.push_back(0);
                                    tOctaCount.push_back(0);
                                }
                                sameCoord = 1;
                            }
    ~ComplexFilterCriteria() {
        alignedQChainTmScores.clear();
        alignedTChainTmScores.clear();
    }

    bool hasTM(float TMThr, int covMode, int filterMode){
        switch (filterMode){
            case LocalParameters::FILTER_MODE_INTERFACE:
            case LocalParameters::FILTER_MODE_LOOSE:
                switch (covMode) {
                    case Parameters::COV_MODE_BIDIRECTIONAL:
                        return ((qTM>= TMThr) && (tTM >= TMThr));
                    case Parameters::COV_MODE_TARGET:   
                        return (tTM >= TMThr);
                    case Parameters::COV_MODE_QUERY:
                        return (qTM >= TMThr);
                    case Parameters::COV_MODE_LENGTH_QUERY :
                    case Parameters::COV_MODE_LENGTH_TARGET :
                    case Parameters::COV_MODE_LENGTH_SHORTER :
                        return true;
                }
            // case LocalParameters::FILTER_MODE_CONFORMATION:
            //TODO 
            //For 1: maybe, check the minimum chain tmscore among all chain-chain(not new) tm
        }
    }

    // bool hasChainNum(int covMode, int filterMode, int qChainNum, int tChainNum ){
    //     switch (filterMode){
    //         case LocalParameters::FILTER_MODE_INTERFACE:
    //             switch (covMode) {
    //                 case Parameters::COV_MODE_BIDIRECTIONAL:
    //                     return (alignedQChainTmScores.size()==qChainNum && qChainNum==tChainNum);
    //                 case Parameters::COV_MODE_TARGET:
    //                     return (alignedTChainTmScores.size()==tChainNum);
    //                 case Parameters::COV_MODE_QUERY:
    //                     return (alignedQChainTmScores.size()==qChainNum);
    //                 case Parameters::COV_MODE_LENGTH_QUERY :
    //                 case Parameters::COV_MODE_LENGTH_TARGET :
    //                 case Parameters::COV_MODE_LENGTH_SHORTER :
    //                     return true;
    //             }
    //         case LocalParameters::FILTER_MODE_CONFORMATION:
    //             switch (covMode) {
    //                 case Parameters::COV_MODE_BIDIRECTIONAL:
    //                     return (qChainNum==tChainNum);
    //                 default:
    //                     return true;
    //             }
    //         case LocalParameters::FILTER_MODE_LOOSE:
    //             return true;

    //     }
    // } 

    bool hasMatchedCoord(int filterMode){
        switch (filterMode) {
            case LocalParameters::FILTER_MODE_INTERFACE:
                return (sameCoord == 1);
            case LocalParameters::FILTER_MODE_CONFORMATION:
            case LocalParameters::FILTER_MODE_LOOSE:
                return true;
        }
    }

    bool hasMatchedOcta(int filterMode, int covMode, std::vector<int> &qOctaCount, std::vector<int> &tOctaCount){
        if (qChainKeys.size()>= 1){
            switch (filterMode) {
                case LocalParameters::FILTER_MODE_INTERFACE:
                    switch (covMode) {
                        case Parameters::COV_MODE_BIDIRECTIONAL:
                            for (size_t i = 0; i < 8; i++) {
                                if (qOctaCount[i] != tOctaCount[i]) {
                                    return false;
                                }
                            }
                        case Parameters::COV_MODE_TARGET:
                            for (size_t i = 0; i < 8; i++) {
                                if (qOctaCount[i] < tOctaCount[i]) {
                                    return false;
                                }
                            }
                            return true;
                        case Parameters::COV_MODE_QUERY:
                            for (size_t i = 0; i < 8; i++) {
                                if (qOctaCount[i] > tOctaCount[i]) {
                                    return false;
                                }
                            }
                            return true; 
                        case Parameters::COV_MODE_LENGTH_QUERY :
                        case Parameters::COV_MODE_LENGTH_TARGET :
                        case Parameters::COV_MODE_LENGTH_SHORTER :
                            return true;
                    }
                case LocalParameters::FILTER_MODE_CONFORMATION:
                case LocalParameters::FILTER_MODE_LOOSE:
                    return true;
            }
        }
        return true;
    }

    bool satisfy(int covMode, int filterMode, float covThr, float TMThr, int qChainNum, int tChainNum ) {
        const bool covOK = Util::hasCoverage(covThr, covMode, qCov, tCov);
        const bool TMOK = hasTM(TMThr, covMode, filterMode);
        // const bool chainNumOK = hasChainNum(covMode, filterMode, qChainNum, tChainNum);
        const bool coordOK = hasMatchedCoord(filterMode);
        // const bool chainTMOK = hasChainTm(chainTMThr, covMode, alignedQChainTmScores, alignedTChainTmScores, qChainNum, tChainNum);
        // return (covOK && TMOK && chainTMOK && numOK);
        const bool octaOK = hasMatchedOcta(filterMode, covMode, qOctaCount, tOctaCount);
        // #pragma omp critical
        // {
        //     if (covOK && TMOK && coordOK){
        //         if(!octaOK){
        //             if (qChainKeys.size()==tChainKeys.size()){
        //                 Debug(Debug::WARNING)<< qChainKeys[0] << "\t"<< tChainKeys[0]<< "\n" ;
        //             }
        //         }
        //     }
        // }
        return (covOK && TMOK && coordOK && octaOK);
        // return(covOK && TMOK && coordOK);
    }

    void update(unsigned int qChainKey, unsigned int tChainKey, unsigned int qTotalAlnLen, unsigned int tTotalAlnLen, double qChainTm, double tChainTm, bool sameCoord, float qAvgCoord[3]) {
        this->tTotalAlnLen += tTotalAlnLen;
        this->qTotalAlnLen += qTotalAlnLen;
        this->alignedQChainTmScores.push_back(qChainTm);
        this->alignedTChainTmScores.push_back(tChainTm);
        this->sameCoord *= sameCoord;
        if (sameCoord){
            refCoord[0] = qAvgCoord[0];
            refCoord[1] = qAvgCoord[1];
            refCoord[2] = qAvgCoord[2];
        }
        auto pos = std::find(qChainKeys.begin(), qChainKeys.end(), qChainKey);
        if (pos != qChainKeys.end()) {
            qChainKeys.erase(pos);
        }
        pos = std::find(tChainKeys.begin(), tChainKeys.end(), tChainKey);
        if (pos != tChainKeys.end()) {
            tChainKeys.erase(pos);
        }
    }

    void calcCov(unsigned int qLen, unsigned int tLen) {
        qCov = static_cast<float>(qTotalAlnLen) / static_cast<float>(qLen);
        tCov = static_cast<float>(tTotalAlnLen) / static_cast<float>(tLen);
    }
};

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

unsigned int cigarToAlignedLength(const std::string &cigar){
    std::string backtrace = Matcher::uncompressAlignment(cigar);
    unsigned int alni = 0;
    for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
        if (backtrace[btPos] == 'M') {
            alni++;
        }
    }
    return alni;
}

void fillMatchedCoord(float * qdata, float * tdata, 
                Coordinates &qm, Coordinates &tm,
                const std::string &cigar, int qStartPos, int tStartPos, int qLen, int tLen) {
    int qi = qStartPos;
    int ti = tStartPos;
    int mi = 0;

    std::string backtrace = Matcher::uncompressAlignment(cigar);
    for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
        if (backtrace[btPos] == 'M') {
            qm.x[mi] = qdata[qi];
            qm.y[mi] = qdata[qLen + qi];
            qm.z[mi] = qdata[2*qLen + qi];
            tm.x[mi] = tdata[ti];
            tm.y[mi] = tdata[tLen + ti];
            tm.z[mi] = tdata[2*tLen + ti];
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
}

double computeChainTmScore(Coordinates &qm, Coordinates &tm, float t[3], float u[3][3], unsigned int alnLen, int tLen) {
    double tmscore = 0;
    float d0 = 1.24*(cbrt(tLen-15)) -1.8;
    float d02 = d0*d0;
    Coordinates tmt(alnLen);
    BasicFunction::do_rotation(tm, tmt, alnLen, t, u);
    for (unsigned int k=0; k<alnLen; k++) {
        double di = BasicFunction::dist(qm.x[k], qm.y[k], qm.z[k], tmt.x[k], tmt.y[k], tmt.z[k]);
        tmscore += 1/(1+di/d02);
    }
    return tmscore;
}

void getRotatedAverageCoord(float * tdata, float tavgCoords[3], int tlen, float t[3], float u[3][3]){
    Coordinates tCoords(tlen);
    for (int i=0; i < tlen ; i++ ){
        tCoords.x[i] = tdata[i];
        tCoords.y[i] = tdata[tlen+i];
        tCoords.z[i] = tdata[2*tlen+i];
    }
    Coordinates Coordsrot(tlen);
    BasicFunction::do_rotation(tCoords, Coordsrot, tlen, t, u);
    float txsum=0.0, tysum=0.0, tzsum=0.0;
    for (int i=0; i< tlen;i++){
        txsum += Coordsrot.x[i];
        tysum += Coordsrot.y[i];
        tzsum += Coordsrot.z[i];
    }
    txsum /= tlen;
    tysum /= tlen;
    tzsum /= tlen;
    tavgCoords[0] = txsum;
    tavgCoords[1] = tysum;
    tavgCoords[2] = tzsum;
}

void getAverageCoord( float * qdata, float qavgCoords[3], int qlen){
    
    float qxsum=0.0, qysum=0.0, qzsum=0.0;
    for (int i=0; i< qlen;i++){
        qxsum += qdata[i];
        qysum += qdata[qlen+i];
        qzsum += qdata[2*qlen+i];
    }
    qxsum /= qlen;
    qysum /= qlen;
    qzsum /= qlen;
    qavgCoords[0] = qxsum;
    qavgCoords[1] = qysum;
    qavgCoords[2] = qzsum;
}

void getVector(float avgCoords[3], float refavgCoords[3]){
    avgCoords[0] -= refavgCoords[0];
    avgCoords[1] -= refavgCoords[1];
    avgCoords[2] -= refavgCoords[2];
}

unsigned int getComplexResidueLength( IndexReader *Dbr, std::vector<unsigned int> &ChainKeys) {
        unsigned int ResidueLen = 0;
        for (auto ChainKey: ChainKeys) {
            size_t id = Dbr->sequenceReader->getId(ChainKey);
            // Not accessible
            if (id == NOT_AVAILABLE_CHAIN_KEY)
                return 0;
            ResidueLen += Dbr->sequenceReader->getSeqLen(id);
        }
        return ResidueLen;
}

static void getlookupInfo(
        const std::string &file,
        std::map<unsigned int, std::string> &complexIdtoName,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::map<unsigned int, std::vector<unsigned int>> &complexIdToChainKeysLookup,
        std::vector<unsigned int> &complexIdVec
) {
    if (file.length() == 0) {
        return;
    }
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];
    int prevComplexId =  -1;
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        auto chainKey = Util::fast_atoi<int>(entry[0]);
        std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
        auto complexId = Util::fast_atoi<int>(entry[2]);
        chainKeyToComplexIdLookup.emplace(chainKey, complexId);
        
        size_t lastUnderscoreIndex = chainName.find_last_of('_');
        std::string complexName = chainName.substr(0, lastUnderscoreIndex);

        if (complexId != prevComplexId) {
            complexIdToChainKeysLookup.emplace(complexId, std::vector<unsigned int>());
            complexIdVec.emplace_back(complexId);
            complexIdtoName.emplace(complexId, complexName);
            prevComplexId = complexId;
        }
        complexIdToChainKeysLookup.at(complexId).emplace_back(chainKey);
        data = Util::skipLine(data);
    }
    lookupDB.close();
}

int filtercomplex(int argc, const char **argv, const Command &command) {
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

    // Debug(Debug::WARNING) << "Monomer will be treated as singleton\nMonomer chain key: \n";
#ifdef OPENMP
localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif
    const bool shouldCompress = (par.compressed == true);
    const int db4Type = Parameters::DBTYPE_CLUSTER_RES;
    
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, shouldCompress, db4Type);
    resultWriter.open();

    const int db5Type = Parameters::DBTYPE_GENERIC_DB;
    DBWriter resultWrite5(par.db5.c_str(), par.db5Index.c_str(), 1, shouldCompress, db5Type);
    resultWrite5.open();

    std::string qLookupFile = par.db1 + ".lookup";
    std::string tLookupFile = par.db2 + ".lookup";
    
    chainKeyToComplexId_t qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;
    complexIdToChainKeys_t qComplexIdToChainKeyMap, tComplexIdToChainKeyMap;
    std::map<unsigned int, std::string> qcomplexIdToName, tcomplexIdToName;
    std::vector<unsigned int> qComplexIdVec, tComplexIdVec;
    getlookupInfo(qLookupFile, qcomplexIdToName,qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
    getlookupInfo(tLookupFile, tcomplexIdToName, tChainKeyToComplexIdMap, tComplexIdToChainKeyMap, tComplexIdVec);
    qChainKeyToComplexIdMap.clear();
    Debug::Progress progress(qComplexIdVec.size());
    std::map<unsigned int, unsigned int> qComplexLength, tComplexLength;
    std::map<unsigned int, std::string> qComplexIdResult;

    for (size_t tComplexIdx = 0; tComplexIdx < tComplexIdVec.size(); tComplexIdx++) {
        unsigned int tComplexId = tComplexIdVec[tComplexIdx];
        std::vector<unsigned int> &tChainKeys = tComplexIdToChainKeyMap.at(tComplexId);
        if (tChainKeys.empty()) {
            continue;
        }
        unsigned int reslen = getComplexResidueLength(tDbr, tChainKeys);
        tComplexLength[tComplexId] =reslen;
    }
    for (size_t qComplexIdx = 0; qComplexIdx < qComplexIdVec.size(); qComplexIdx++) {
        unsigned int qComplexId = qComplexIdVec[qComplexIdx];
        std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap.at(qComplexId);
        if (qChainKeys.empty()) {
            continue;
        }
        unsigned int reslen = getComplexResidueLength(qDbr, qChainKeys);
        qComplexLength[qComplexId] = reslen;
    }

    
    
#pragma omp parallel num_threads(localThreads) 
    {   
        resultToWrite_t result5;
        char buffer[32];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string result;
        std::map<unsigned int, std::string> tmpDBKEYut;
        std::map<unsigned int, ComplexFilterCriteria> localComplexMap;
        std::vector<unsigned int> assIdsToDelete;
        std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId; // cmplId : [assId, alnSum]
        std::vector<unsigned int> selectedAssIDs;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        
        Matcher::result_t res;
#pragma omp for schedule(dynamic, 1)    
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            progress.updateProgress();
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
            std::vector<unsigned int> qChainKeys = qComplexIdToChainKeyMap.at(qComplexId);

            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainAlnId = alnDbr.getId(qChainKey);
                unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                //handling monomer as singleton
                if (qChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                    char *outpos = Itoa::u32toa_sse2(qComplexId, buffer);
                    result.append(buffer, (outpos - buffer - 1));
                    result.push_back('\n');
                    result5.append(qcomplexIdToName.at(qComplexId) + "\t" + tcomplexIdToName.at(qComplexId) + "\t1.000000\t1.000000\t1.000000\t1.000000\n");
                    break;
                }
                
                char *data = alnDbr.getData(qChainAlnId, thread_idx);
                while (*data != '\0' ) {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
                    size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
                    float* qdata = qcoords.read(qcadata, res.qLen, qCaLength);
                
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }

                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainKey= res.dbKey;
                    unsigned int tChainAlnId = alnDbr.getId(tChainKey);
                    //if target is monomer, break to be singleton
                    if (tChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                        break;
                    }
                    unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                    unsigned int tComplexId = tChainKeyToComplexIdMap.at(tChainKey);
                    std::vector<unsigned int> tChainKeys = tComplexIdToChainKeyMap.at(tComplexId);
                    
                    float u[3][3];
                    float t[3];
                    fillUArr(retComplex.uString, u);
                    fillTArr(retComplex.tString, t);
                    tmpDBKEYut[assId]=retComplex.uString+","+retComplex.tString;
                    char *tcadata = tStructDbr->getData(tChainDbId, thread_idx);
                    size_t tCaLength = tStructDbr->getEntryLen(tChainDbId);
                    float* tdata = tcoords.read(tcadata, res.dbLen, tCaLength);
                    unsigned int alnLen = cigarToAlignedLength(res.backtrace);
                    Coordinates qm(alnLen), tm(alnLen);
                    //FIXME: if new chainTM not required, erase those part
                    fillMatchedCoord(qdata, tdata, qm, tm, res.backtrace, res.qStartPos, res.dbStartPos, res.qLen, res.dbLen);
                    double chainTm = computeChainTmScore(qm, tm, t, u, alnLen, res.dbLen);
                    double qChainTm = chainTm / res.qLen;
                    double tChainTm = chainTm/ res.dbLen;
                    unsigned int qtotalaln = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    unsigned int ttotalaln = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                    
                    float qAvgCoord[3], tAvgCoord[3];
                    getRotatedAverageCoord(tdata, tAvgCoord, res.dbLen, t, u);
                    getAverageCoord(qdata, qAvgCoord, res.qLen) ;
                    bool coordSame = isSameCoord(qAvgCoord, tAvgCoord);
                    if (localComplexMap.find(assId) == localComplexMap.end()) {
                        ComplexFilterCriteria cmplfiltcrit = ComplexFilterCriteria(tChainKey, qChainKeys, tChainKeys, retComplex.qTmScore, retComplex.tTmScore, t, u);
                        localComplexMap[assId] = cmplfiltcrit;
                        localComplexMap.at(assId).update(qChainKey, tChainKey, qtotalaln, ttotalaln, qChainTm, tChainTm, coordSame, qAvgCoord);
                    } else {
                        localComplexMap.at(assId).update(qChainKey, tChainKey, qtotalaln, ttotalaln, qChainTm, tChainTm, coordSame, qAvgCoord);
                    }
                    
                } // while end
            }
            for (auto& assId_res : localComplexMap){
                unsigned int tComplexId = tChainKeyToComplexIdMap.at(assId_res.second.dbKey);
                assId_res.second.calcCov(qComplexLength.at(qComplexId), tComplexLength.at(tComplexId));
                std::vector<unsigned int> tChainKeys = tComplexIdToChainKeyMap.at(tComplexId);
                if (assId_res.second.qChainKeys.size()!= 0){
                    std::vector<int> qOctaCount(8,0);
                    std::vector<int> tOctaCount(8,0);
                    for (unsigned int qChainKey : assId_res.second.qChainKeys){
                        unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                        char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
                        size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
                        int qSeqLength = qDbr->sequenceReader->getSeqLen(qChainDbId);
                        float* qdata = qcoords.read(qcadata, qSeqLength, qCaLength);
                        float qAvgCoord[3];
                        getAverageCoord(qdata, qAvgCoord, qSeqLength);
                        int tocheck = determineOctant(qAvgCoord);
                        getVector(qAvgCoord, assId_res.second.refCoord);
                        int qOcta = determineOctant(qAvgCoord);
                        qOctaCount[qOcta]++;
                    }
                    for (unsigned int tChainKey : assId_res.second.tChainKeys){
                        unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                        char *tcadata = tStructDbr->getData(tChainDbId, thread_idx);
                        size_t tCaLength = tStructDbr->getEntryLen(tChainDbId);
                        int tSeqLength = tDbr->sequenceReader->getSeqLen(tChainDbId);
                        float* tdata = tcoords.read(tcadata, tSeqLength, tCaLength);
                        float tAvgCoord[3];
                        getRotatedAverageCoord(tdata, tAvgCoord, tSeqLength, assId_res.second.t, assId_res.second.u);
                        getVector(tAvgCoord, assId_res.second.refCoord);
                        int tOcta = determineOctant(tAvgCoord);
                        tOctaCount[tOcta]++;
                    }
                    for (int i=0;i<8;i++){
                        assId_res.second.qOctaCount[i] = qOctaCount[i];
                        assId_res.second.tOctaCount[i] = tOctaCount[i];
                    }
                }

                if (!(assId_res.second.satisfy(par.covMode, par.filterMode, par.covThr, par.filtComplexTmThr, qChainKeys.size(), tChainKeys.size()))){
                    assIdsToDelete.push_back(assId_res.first);
                }
            }

            for (const auto& key : assIdsToDelete) {
                localComplexMap.erase(key);
            }
            
            for (const auto& assId_res : localComplexMap){
                unsigned int tComplexId = tChainKeyToComplexIdMap.at(assId_res.second.dbKey);
                unsigned int alnlen = adjustAlnLen(assId_res.second.qTotalAlnLen, assId_res.second.tTotalAlnLen, par.covMode);
                if (cmplIdToBestAssId.find(tComplexId) == cmplIdToBestAssId.end()){
                    cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                }
                else {
                    if (alnlen > cmplIdToBestAssId.at(tComplexId)[1]){
                        cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    }
                }
            }

            for (const auto& pair : cmplIdToBestAssId){
                selectedAssIDs.push_back(pair.second[0]);
            }

            for (unsigned int assIdidx = 0; assIdidx < selectedAssIDs.size(); assIdidx++){
                unsigned int assId = selectedAssIDs[assIdidx];
                unsigned int tComplexId = tChainKeyToComplexIdMap.at(localComplexMap.at(assId).dbKey);
                char *outpos = Itoa::u32toa_sse2(tComplexId, buffer);
                result.append(buffer, (outpos - buffer - 1));
                result.push_back('\n');
                result5.append(qcomplexIdToName.at(qComplexId) + "\t" + tcomplexIdToName.at(tComplexId) + "\t" + std::to_string(localComplexMap.at(assId).qCov) + "\t" + std::to_string(localComplexMap.at(assId).tCov) + "\t"+ std::to_string(localComplexMap.at(assId).qTM)+"\t"+ std::to_string(localComplexMap.at(assId).tTM)+ "\n");
            }
            #pragma omp critical
            {
                qComplexIdResult[qComplexId]= result;
            }
            result.clear();
            localComplexMap.clear();
            tmpDBKEYut.clear();
            assIdsToDelete.clear();
            cmplIdToBestAssId.clear();
            selectedAssIDs.clear();
        } // for end
        #pragma omp critical
        {
            resultWrite5.writeData(result5.c_str(), result5.length(), 0);
        }
        result5.clear();
    } // MP end
    for (auto &pair : qComplexIdResult){
        resultWriter.writeData(pair.second.c_str(), pair.second.length(), pair.first);
    }
    
    resultWriter.close(true);
    resultWrite5.close(par.dbOut == false);
    qStructDbr.close();
    alnDbr.close();
    delete qDbr;
    if (sameDB == false) {
        delete tDbr;
        delete tStructDbr;
    }
    qChainKeyToComplexIdMap.clear();
    tChainKeyToComplexIdMap.clear();
    qComplexIdToChainKeyMap.clear();
    tComplexIdToChainKeyMap.clear();
    qcomplexIdToName.clear();
    tcomplexIdToName.clear();
    qComplexIdVec.clear();
    tComplexIdVec.clear();
    qComplexLength.clear();
    tComplexLength.clear();
    
    return EXIT_SUCCESS;
}
