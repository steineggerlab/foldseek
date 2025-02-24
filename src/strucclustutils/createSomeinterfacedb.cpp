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

//one dimer db as an input, one interface db as an output
int createSomeinterfacedb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<unsigned int> qDbr((par.db1).c_str(), (par.db1 + ".index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qDbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qStructDbr.open(DBReader<unsigned int>::NOSORT);

    std::string qLookupFile = par.db1 + ".lookup";
    std::vector<unsigned int> qComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    qChainKeyToComplexIdMap.clear();
    
    DBWriter ssdbw((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    ssdbw.open();
    DBWriter cadbw((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    cadbw.open();
    DBWriter aadbw((par.db2).c_str(), (par.db2 + ".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aadbw.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
    #endif
        std::vector<int8_t> camol1, camol2;
        std::vector<float> ca1, ca2;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
#pragma omp for schedule(static)
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            if (qChainKeys.size() != 2) {
                continue;
            }
            unsigned int qChainIdx = 0;
            unsigned int qChainKey = qChainKeys[qChainIdx];
            unsigned int qChainDbId = qDbr.getId(qChainKey);
            char *qaaadata = qDbr.getData(qChainDbId, thread_idx);
            char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
            size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
            size_t qChainLen = qDbr.getSeqLen(qChainDbId);
            float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);

            unsigned int tChainIdx = 1;
            unsigned int tChainKey = qChainKeys[tChainIdx];
            unsigned int tChainDbId = qDbr.getId(tChainKey);
            char *taaadata = qDbr.getData(tChainDbId, thread_idx);
            char *tcadata = qStructDbr.getData(tChainDbId, thread_idx);
            size_t tCaLength = qStructDbr.getEntryLen(tChainDbId);
            size_t tChainLen = qDbr.getSeqLen(tChainDbId);
            float* tdata = qcoords.read(tcadata, tChainLen, tCaLength);
            
            float distanceThreshold = 10;
            std::vector<size_t> resIdx1, resIdx2;
            const double squareThreshold = distanceThreshold * distanceThreshold;
            PulchraWrapper pulchra;
            std::vector<Vec3> caB, nB, cB;
            std::vector<char> amiB;
            amiB.resize(qChainLen + tChainLen);
            caB.resize(qChainLen + tChainLen);
            nB.resize(qChainLen + tChainLen);
            cB.resize(qChainLen + tChainLen);
            for (size_t i = 0; i < qChainLen; i++) {
                caB[i] = Vec3(qdata[i], qdata[i + qChainLen], qdata[i + qChainLen * 2]);
                nB[i] = Vec3(0,0,0);
                cB[i] = Vec3(0,0,0);
                amiB[i] = qaaadata[i];
            }
            for (size_t i = 0; i < tChainLen; i++) {
                caB[qChainLen + i] = Vec3(tdata[i], tdata[i + tChainLen], tdata[i+ tChainLen * 2]);
                nB[qChainLen + i] = Vec3(0,0,0);
                cB[qChainLen + i] = Vec3(0,0,0);
                amiB[qChainLen + i] = taaadata[i];
            }
            pulchra.rebuildBackbone(&caB[0], &nB[0], &cB[0], &amiB[0], qChainLen + tChainLen);
            findInterface(resIdx1, squareThreshold, qdata, tdata, qChainLen, tChainLen);
            findInterface(resIdx2, squareThreshold, tdata, qdata, tChainLen, qChainLen);
            if (resIdx1.size() >= 4 && resIdx2.size() >= 4) {
                // std::cout<<qChainKey<<"\t"<<tChainKey<<std::endl;
                StructureTo3Di structureTo3Di;
                std::vector<Vec3> ca, n, c, cb;
                std::vector<char> ami;
                std::vector<char> alphabet3di1, alphabet3di2;
                std::vector<char> alphabetAA1, alphabetAA2;
                SubstitutionMatrix mat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
                ami.resize(resIdx1.size() + resIdx2.size());
                ca.resize(resIdx1.size() + resIdx2.size());
                n.resize(resIdx1.size() + resIdx2.size());
                c.resize(resIdx1.size() + resIdx2.size());
                cb.resize(resIdx1.size() + resIdx2.size());
                for (size_t i = 0; i < resIdx1.size(); i++) {
                    ca[i] = caB[resIdx1[i]];
                    n[i] = nB[resIdx1[i]];
                    c[i] = cB[resIdx1[i]];
                    cb[i] = Vec3(0,0,0);
                    ami[i] = amiB[resIdx1[i]];
                }
                for (size_t i = 0; i < resIdx2.size(); i++) {
                    ca[resIdx1.size() + i] = caB[qChainLen + resIdx2[i]];
                    n[resIdx1.size() + i] = nB[qChainLen + resIdx2[i]];
                    c[resIdx1.size() + i] = cB[qChainLen + resIdx2[i]];
                    cb[resIdx1.size() + i] = Vec3(0,0,0);
                    ami[resIdx1.size() + i] = amiB[qChainLen + resIdx2[i]];
                }

                char *states = structureTo3Di.structure2states(&ca[0],
                                                                    &n[0],
                                                                    &c[0],
                                                                    &cb[0],
                                                                    resIdx1.size() + resIdx2.size());
                alphabet3di1.resize(resIdx1.size() + 1);
                alphabetAA1.resize(resIdx1.size() + 1);
                alphabet3di2.resize(resIdx2.size() + 1);                               
                alphabetAA2.resize(resIdx2.size() + 1);                             
                ca1.resize(3 * resIdx1.size());                          
                ca2.resize(3 * resIdx2.size());
                camol1.resize((resIdx1.size() - 1) * 3 * sizeof(int16_t) + 3 * sizeof(float));
                camol2.resize((resIdx2.size() - 1) * 3 * sizeof(int16_t) + 3 * sizeof(float));
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
                // ca1[resIdx1.size()] = '\n';
                for (size_t pos = 0; pos < resIdx2.size(); pos++) {
                    alphabet3di2[pos] = mat.num2aa[static_cast<int>(states[resIdx1.size() + pos])];
                    alphabetAA2[pos] = ami[resIdx1.size() + pos];
                    ca2[pos] = ca[resIdx1.size() + pos].x;
                    ca2[pos + resIdx2.size()] = ca[resIdx1.size() + pos].y;
                    ca2[pos + 2 * resIdx2.size()] = ca[resIdx1.size() + pos].z;
                }
                alphabet3di2[resIdx2.size()] = '\n';
                alphabetAA2[resIdx2.size()] = '\n';
                // ca2[resIdx2.size()] = '\n';
                ssdbw.writeData(alphabet3di1.data(), alphabet3di1.size(), qChainKey, thread_idx);
                ssdbw.writeData(alphabet3di2.data(), alphabet3di2.size(), tChainKey, thread_idx);
                aadbw.writeData(alphabetAA1.data(), alphabetAA1.size(), qChainKey, thread_idx);
                aadbw.writeData(alphabetAA2.data(), alphabetAA2.size(), tChainKey, thread_idx);
                char *data1 = reinterpret_cast<char*>(ca1.data());
                char *data2 = reinterpret_cast<char*>(ca2.data());
                if (!Coordinate16::convertToDiff16(resIdx1.size(), (float*)(data1), camol1f16, 1)
                        && !Coordinate16::convertToDiff16(resIdx1.size(), (float*)(data1) + resIdx1.size(), camol1f16 + resIdx1.size(), 1)
                        && !Coordinate16::convertToDiff16(resIdx1.size(), (float*)(data1) + 2 * resIdx1.size(), camol1f16 + 2 * resIdx1.size(), 1)) {
                    cadbw.writeData((const char*)camol1.data(), (resIdx1.size() - 1) * 3 * sizeof(uint16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), qChainKey, thread_idx);
                } else {
                    cadbw.writeData(data1, resIdx1.size() * 3 * sizeof(float), qChainKey, thread_idx);
                }
                if (!Coordinate16::convertToDiff16(resIdx2.size(), (float*)(data2), camol2f16, 1)
                        && !Coordinate16::convertToDiff16(resIdx2.size(), (float*)(data2) + resIdx2.size(), camol2f16 + resIdx2.size(), 1)
                        && !Coordinate16::convertToDiff16(resIdx2.size(), (float*)(data2) + 2 * resIdx2.size(), camol2f16 + 2 * resIdx2.size(), 1)) {
                    cadbw.writeData((const char*)camol2.data(), (resIdx2.size() - 1) * 3 * sizeof(uint16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), tChainKey, thread_idx);
                } else {
                    cadbw.writeData(data2, resIdx2.size() * 3 * sizeof(float), tChainKey, thread_idx);
                }
            }
            ca1.clear();
            ca2.clear();
            camol1.clear();
            camol2.clear();
        }
    }
    ssdbw.close(true);
    cadbw.close(true);
    aadbw.close(true);
    qStructDbr.close();
    qDbr.close();

    return EXIT_SUCCESS;
}
