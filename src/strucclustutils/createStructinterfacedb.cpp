#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "tmalign/basic_fun.h"
#include "MultimerUtil.h"
#include "structureto3di.h"
#include "SubstitutionMatrix.h"
#include "GemmiWrapper.h"
#include "PulchraWrapper.h"
#include "LDDT.h"
#include "itoa.h"
#include "MathUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

// read one residue index out of the packed _id entry, without copying the whole chain
static inline unsigned int getResId(const char *idData, size_t resIdx) {
    unsigned int resId;
    memcpy(&resId, idData + resIdx * sizeof(unsigned int), sizeof(unsigned int));
    return resId;
}

// pulchra needs two residues of context on each side of a residue it reconstructs, and
// prepare_rbins extrapolates a window's ends from its first/last five residues
const size_t PULCHRA_HALO = 2;
const size_t PULCHRA_MIN_WINDOW = 5;

//one dimer db as an input, one interface db as an output
int createStructinterfacedb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    bool saveResIndex = true;
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<unsigned int> qDbr((par.db1).c_str(), (par.db1 + ".index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qDbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qStructDbr.open(DBReader<unsigned int>::NOSORT);
    std::string idDbFileName = par.db1;
    idDbFileName.append("_id");
    const char *idDbChar = idDbFileName.c_str();
    std::vector<std::string> dataFileNames = FileUtil::findDatafiles(idDbChar);
    if (dataFileNames.empty()) {
        saveResIndex = false;
    }
    dataFileNames.clear();
    DBReader<unsigned int>* qIdDbr = NULL;
    if (saveResIndex) {
        qIdDbr = new DBReader<unsigned int>((par.db1 + "_id").c_str(), (par.db1 + "_id.index").c_str(),
                            par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        qIdDbr->open(DBReader<unsigned int>::NOSORT);
    }

    std::string qLookupFile = par.db1 + ".lookup";
    std::vector<unsigned int> qComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    qChainKeyToComplexIdMap.clear();

    // flatten the dimers into a contiguous list so that the parallel loop needs no map lookup,
    // is not unbalanced by skipped complexes, and the lookup maps can be released up front
    std::vector<std::pair<unsigned int, unsigned int>> dimerChainKeys;
    dimerChainKeys.reserve(qComplexIndices.size());
    for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
        const std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexIndices[qCompIdx]);
        if (qChainKeys.size() != 2) {
            continue;
        }
        dimerChainKeys.emplace_back(qChainKeys[0], qChainKeys[1]);
    }
    qComplexIdToChainKeysMap.clear();
    qComplexIndices.clear();
    qComplexIndices.shrink_to_fit();

    int aadbType = Parameters::DBTYPE_AMINO_ACIDS;
    int cadbType = LocalParameters::DBTYPE_CA_ALPHA;

    aadbType = DBReader<unsigned int>::setExtendedDbtype(aadbType, LocalParameters::DBTYPE_EXTENDED_INTERFACE);
    cadbType = DBReader<unsigned int>::setExtendedDbtype(cadbType, LocalParameters::DBTYPE_EXTENDED_INTERFACE);

    DBWriter ssdbw((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, aadbType);
    ssdbw.open();
    DBWriter cadbw((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, cadbType);
    cadbw.open();
    DBWriter aadbw((par.db2).c_str(), (par.db2 + ".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, aadbType);
    aadbw.open();
    DBWriter* idbw = NULL;

    if (saveResIndex) {
        int iddbType = Parameters::DBTYPE_GENERIC_DB;
        iddbType = DBReader<unsigned int>::setExtendedDbtype(iddbType, LocalParameters::DBTYPE_EXTENDED_INTERFACE);
        idbw = new DBWriter((par.db2 + "_id").c_str(), (par.db2 + "_id.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, iddbType);
        idbw->open();
    }
    float distanceThreshold = par.distanceThreshold;
    unsigned int minimumResidue = par.minResidueNum;
    const float squareThreshold = distanceThreshold * distanceThreshold;

    // only num2aa is needed, parsing the matrix once is enough for all threads
    SubstitutionMatrix mat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    Debug::Progress progress(dimerChainKeys.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
    #endif
        // these are expensive to construct (pulchra allocates per-residue scratch,
        // StructureTo3Di loads the kerasify model), so keep one per thread
        StructureTo3Di structureTo3Di;
        PulchraWrapper pulchra;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        std::vector<int8_t> camol1, camol2;
        std::vector<float> ca1, ca2;
        std::vector<size_t> resIdx1, resIdx2;
        std::vector<Vec3> caB, nB, cB;
        std::vector<char> amiB;
        std::vector<Vec3> ca, n, c, cb;
        std::vector<char> ami;
        std::vector<char> alphabet3di1, alphabet3di2;
        std::vector<char> alphabetAA1, alphabetAA2;
        std::vector<unsigned int> resId1, resId2;
        std::vector<std::pair<size_t, size_t>> pulchraWindows;
        std::vector<unsigned char> qFlag, tFlag;
        // dynamic: per-dimer cost is O(qLen*tLen) and varies by orders of magnitude.
        // note this makes the physical entry order in the data file depend on thread
        // timing; entry content and the key-sorted index are unaffected
#pragma omp for schedule(dynamic, 1)
        for (size_t dimerIdx = 0; dimerIdx < dimerChainKeys.size(); dimerIdx++) {
            progress.updateProgress();
            unsigned int qChainKey = dimerChainKeys[dimerIdx].first;
            unsigned int qChainDbId = qDbr.getId(qChainKey);
            char *qaaadata = qDbr.getData(qChainDbId, thread_idx);
            char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
            size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
            size_t qChainLen = qDbr.getSeqLen(qChainDbId);
            float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);

            unsigned int tChainKey = dimerChainKeys[dimerIdx].second;
            unsigned int tChainDbId = qDbr.getId(tChainKey);
            char *taaadata = qDbr.getData(tChainDbId, thread_idx);
            char *tcadata = qStructDbr.getData(tChainDbId, thread_idx);
            size_t tCaLength = qStructDbr.getEntryLen(tChainDbId);
            size_t tChainLen = qDbr.getSeqLen(tChainDbId);
            float* tdata = tcoords.read(tcadata, tChainLen, tCaLength);

            // the interface only depends on the C-alpha coordinates, so decide whether this
            // dimer is kept before paying for the backbone reconstruction
            findInterface(resIdx1, resIdx2, squareThreshold, qdata, tdata, qChainLen, tChainLen, qFlag, tFlag);
            if (resIdx1.size() < minimumResidue || resIdx2.size() < minimumResidue) {
                continue;
            }

            const size_t chainLenSum = qChainLen + tChainLen;
            amiB.resize(chainLenSum);
            caB.resize(chainLenSum);
            nB.resize(chainLenSum);
            cB.resize(chainLenSum);
            for (size_t i = 0; i < qChainLen; i++) {
                caB[i] = Vec3(qdata[i], qdata[i + qChainLen], qdata[i + qChainLen * 2]);
                nB[i] = Vec3(NAN,NAN,NAN);
                cB[i] = Vec3(NAN,NAN,NAN);
                amiB[i] = qaaadata[i];
            }
            for (size_t i = 0; i < tChainLen; i++) {
                caB[qChainLen + i] = Vec3(tdata[i], tdata[i + tChainLen], tdata[i+ tChainLen * 2]);
                nB[qChainLen + i] = Vec3(NAN,NAN,NAN);
                cB[qChainLen + i] = Vec3(NAN,NAN,NAN);
                amiB[qChainLen + i] = taaadata[i];
            }
            // pulchra reconstructs residue j purely from ca[j-2..j+2]: prepare_rbins only
            // takes the 1-3/1-4 distances of that window and the main loop superimposes a
            // template picked from the fixed nco_stat table, so there is no spatial neighbour
            // lookup over the residues. Rebuilding merged windows around the interface
            // residues therefore yields the same N/C as rebuilding both chains in full.
            // Windows live in the concatenated index space, so behaviour at the junction of
            // the two chains is unchanged.
            pulchraWindows.clear();
            for (int ch = 0; ch < 2; ch++) {
                const std::vector<size_t> &res = (ch == 0) ? resIdx1 : resIdx2;
                const size_t base = (ch == 0) ? 0 : qChainLen;
                for (size_t i = 0; i < res.size(); i++) {
                    const size_t idx = base + res[i];
                    const size_t begin = (idx > PULCHRA_HALO) ? idx - PULCHRA_HALO : 0;
                    size_t end = idx + PULCHRA_HALO + 1;
                    if (end > chainLenSum) {
                        end = chainLenSum;
                    }
                    if (pulchraWindows.empty() == false && begin <= pulchraWindows.back().second) {
                        if (end > pulchraWindows.back().second) {
                            pulchraWindows.back().second = end;
                        }
                    } else {
                        pulchraWindows.emplace_back(begin, end);
                    }
                }
            }
            for (size_t w = 0; w < pulchraWindows.size(); w++) {
                size_t begin = pulchraWindows[w].first;
                size_t end = pulchraWindows[w].second;
                // prepare_rbins extrapolates the window ends from five residues, so never
                // hand it a shorter window than that
                while (end - begin < PULCHRA_MIN_WINDOW && end < chainLenSum) {
                    end++;
                }
                while (end - begin < PULCHRA_MIN_WINDOW && begin > 0) {
                    begin--;
                }
                pulchra.rebuildBackbone(&caB[begin], &nB[begin], &cB[begin], &amiB[begin], end - begin);
            }

            const size_t interfaceLen = resIdx1.size() + resIdx2.size();
            ami.resize(interfaceLen);
            ca.resize(interfaceLen);
            n.resize(interfaceLen);
            c.resize(interfaceLen);
            cb.resize(interfaceLen);
            if (saveResIndex) {
                resId1.resize(resIdx1.size());
                resId2.resize(resIdx2.size());
            }
            const char *qiddata = saveResIndex ? qIdDbr->getData(qChainDbId, thread_idx) : NULL;
            const char *tiddata = saveResIndex ? qIdDbr->getData(tChainDbId, thread_idx) : NULL;
            for (size_t i = 0; i < resIdx1.size(); i++) {
                ca[i] = caB[resIdx1[i]];
                n[i] = nB[resIdx1[i]];
                c[i] = cB[resIdx1[i]];
                cb[i] = Vec3(NAN,NAN,NAN);
                ami[i] = amiB[resIdx1[i]];
                if (saveResIndex) {
                    resId1[i] = getResId(qiddata, resIdx1[i]);
                }
            }
            for (size_t i = 0; i < resIdx2.size(); i++) {
                ca[resIdx1.size() + i] = caB[qChainLen + resIdx2[i]];
                n[resIdx1.size() + i] = nB[qChainLen + resIdx2[i]];
                c[resIdx1.size() + i] = cB[qChainLen + resIdx2[i]];
                cb[resIdx1.size() + i] = Vec3(NAN,NAN,NAN);
                ami[resIdx1.size() + i] = amiB[qChainLen + resIdx2[i]];
                if (saveResIndex) {
                    resId2[i] = getResId(tiddata, resIdx2[i]);
                }
            }

            char *states = structureTo3Di.structure2states(&ca[0],
                                                                &n[0],
                                                                &c[0],
                                                                &cb[0],
                                                                interfaceLen);
            alphabet3di1.resize(resIdx1.size() + 1);
            alphabetAA1.resize(resIdx1.size() + 1);
            alphabet3di2.resize(resIdx2.size() + 1);
            alphabetAA2.resize(resIdx2.size() + 1);
            ca1.resize(3 * resIdx1.size());
            ca2.resize(3 * resIdx2.size());
            // assign, not resize: convertToDiff16 leaves the trailing padding byte untouched
            // and it is part of the written entry, so it has to be zeroed like everywhere else
            camol1.assign((resIdx1.size() - 1) * 3 * sizeof(int16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), 0);
            camol2.assign((resIdx2.size() - 1) * 3 * sizeof(int16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), 0);
            int16_t* camol1f16 = reinterpret_cast<int16_t*>(camol1.data());
            int16_t* camol2f16 = reinterpret_cast<int16_t*>(camol2.data());
            for (size_t pos = 0; pos < resIdx1.size(); pos++) {
                alphabet3di1[pos] = mat.num2aa[static_cast<int>(states[pos])];
                alphabetAA1[pos] = ami[pos];
                ca1[pos] = ca[pos].x;
                ca1[pos + resIdx1.size()] = ca[pos].y;
                ca1[pos + 2 * resIdx1.size()] = ca[pos].z;
            }
            alphabet3di1[resIdx1.size()] = '\n';
            alphabetAA1[resIdx1.size()] = '\n';
            for (size_t pos = 0; pos < resIdx2.size(); pos++) {
                alphabet3di2[pos] = mat.num2aa[static_cast<int>(states[resIdx1.size() + pos])];
                alphabetAA2[pos] = ami[resIdx1.size() + pos];
                ca2[pos] = ca[resIdx1.size() + pos].x;
                ca2[pos + resIdx2.size()] = ca[resIdx1.size() + pos].y;
                ca2[pos + 2 * resIdx2.size()] = ca[resIdx1.size() + pos].z;
            }
            alphabet3di2[resIdx2.size()] = '\n';
            alphabetAA2[resIdx2.size()] = '\n';
            ssdbw.writeData(alphabet3di1.data(), alphabet3di1.size(), qChainKey, thread_idx);
            ssdbw.writeData(alphabet3di2.data(), alphabet3di2.size(), tChainKey, thread_idx);
            aadbw.writeData(alphabetAA1.data(), alphabetAA1.size(), qChainKey, thread_idx);
            aadbw.writeData(alphabetAA2.data(), alphabetAA2.size(), tChainKey, thread_idx);
            if (saveResIndex) {
                idbw->writeData((const char*)resId1.data(), resIdx1.size() * sizeof(unsigned int), qChainKey, thread_idx);
                idbw->writeData((const char*)resId2.data(), resIdx2.size() * sizeof(unsigned int), tChainKey, thread_idx);
            }
            char *data1 = reinterpret_cast<char*>(ca1.data());
            char *data2 = reinterpret_cast<char*>(ca2.data());
            if (!Coordinate16::convertToDiff16(resIdx1.size(), (float*)(data1), camol1f16, 1)
                    && !Coordinate16::convertToDiff16(resIdx1.size(), (float*)(data1) + resIdx1.size(), camol1f16 + 1 * (resIdx1.size() + 1), 1)
                    && !Coordinate16::convertToDiff16(resIdx1.size(), (float*)(data1) + 2 * resIdx1.size(), camol1f16 + 2 * (resIdx1.size() + 1), 1)) {
                cadbw.writeData((const char*)camol1.data(), (resIdx1.size() - 1) * 3 * sizeof(uint16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), qChainKey, thread_idx);
            } else {
                cadbw.writeData(data1, resIdx1.size() * 3 * sizeof(float), qChainKey, thread_idx);
            }
            if (!Coordinate16::convertToDiff16(resIdx2.size(), (float*)(data2), camol2f16, 1)
                    && !Coordinate16::convertToDiff16(resIdx2.size(), (float*)(data2) + resIdx2.size(), camol2f16 + 1 * (resIdx2.size() + 1), 1)
                    && !Coordinate16::convertToDiff16(resIdx2.size(), (float*)(data2) + 2 * resIdx2.size(), camol2f16 + 2 * (resIdx2.size() + 1), 1)) {
                cadbw.writeData((const char*)camol2.data(), (resIdx2.size() - 1) * 3 * sizeof(uint16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), tChainKey, thread_idx);
            } else {
                cadbw.writeData(data2, resIdx2.size() * 3 * sizeof(float), tChainKey, thread_idx);
            }
        }
    }
    ssdbw.close(true);
    cadbw.close(true);
    aadbw.close(true);
    qStructDbr.close();
    qDbr.close();
    if (saveResIndex) {
        idbw->close(true);
        delete idbw;
        qIdDbr->close();
        delete qIdDbr;
    }
    return EXIT_SUCCESS;
}
