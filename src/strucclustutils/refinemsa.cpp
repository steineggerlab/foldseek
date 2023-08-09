#include <fstream>
#include <iostream>
#include <vector>
#include "DBReader.h"
#include "Sequence.h"
#include "kseq.h"
#include "KSeqBufferReader.h"
#include "KSeqWrapper.h"
#include "Matcher.h"
#include "structuremsa.h"
#include "msa2lddt.h"
#include "LocalParameters.h"
#include "IndexReader.h"
#include "DBWriter.h"

void deleteGapCols(std::vector<std::string> &sequences) {
    int j = sequences[0].length() - 1;
    while (j >= 0) {
        bool delCol = true;
        for (size_t k = 0; k < sequences.size(); k++) {
            if (sequences[k][j] != '-') {
                delCol = false;
                break;
            }
        }            
        if (delCol) {
            for (size_t k = 0; k < sequences.size(); k++)
                sequences[k].erase(j, 1);
        }
        j--;
    }
}

void buildSubMSA(std::vector<std::string> &headers, std::vector<std::string> &sequences, std::string &subMSA) {
    for (size_t j = 0; j < headers.size(); j++) {
        subMSA.append(1, '>');
        subMSA.append(headers[j]);
        subMSA.append(1, '\n');
        subMSA.append(sequences[j]);
        subMSA.append(1, '\n');
    }
}

void makeSubMSA(std::string msa, std::string &subMSA1, std::string &subMSA2, std::vector<bool> &group) {
    std::vector<std::string> subHeaders1;
    std::vector<std::string> subHeaders2;
    std::vector<std::string> subSeqs1;
    std::vector<std::string> subSeqs2;
    int i = 0;
    KSeqBuffer *kseq = new KSeqBuffer(msa.c_str(), msa.length());
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;
        std::string seq = std::string(entry.sequence.s, entry.sequence.l);
        if (group[i]) {
            subHeaders1.push_back(entry.name.s);
            subSeqs1.push_back(seq);
        } else {
            subHeaders2.push_back(entry.name.s);
            subSeqs2.push_back(seq);
        }
        i++;
    }
    delete kseq;
    deleteGapCols(subSeqs1);
    deleteGapCols(subSeqs2);
    buildSubMSA(subHeaders1, subSeqs1, subMSA1);
    buildSubMSA(subHeaders2, subSeqs2, subMSA2);
}

/**
 * @brief Converts an amino acid MSA to another alphabet using a Foldseek DB.
 * 
 * @param seqDbr3Di 3Di DBReader
 * @param msaAa Amino acid MSA
 * @return std::string 3Di MSA
 */
std::string convertMSA(DBReader<unsigned int> &seqDbrAA, DBReader<unsigned int> &seqDbr, std::string msa) {
    std::string msaNew;
    
    KSeqBuffer *kseq = new KSeqBuffer(msa.c_str(), msa.length());
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;

        // Get corresponding 3Di sequence
        size_t idx = seqDbrAA.getLookupIdByAccession(entry.name.s);
        unsigned int key = seqDbrAA.getLookupKey(idx);
        size_t seqId = seqDbrAA.getId(key);
        std::string seq = std::string(seqDbr.getData(seqId, 0), seqDbr.getSeqLen(seqId));

        // Convert AA to 3Di
        std::string name = std::string(entry.name.s, entry.name.l);
        std::string seqAa = std::string(entry.sequence.s, entry.sequence.l);
        msaNew.push_back('>');
        msaNew.append(name);
        msaNew.push_back('\n');
        int i = 0;
        for (char c : seqAa) {
            if (c == '-') {
                msaNew.push_back(c);        
            } else {
                msaNew.push_back(seq[i]);
                i++;
            }
        }
        msaNew.push_back('\n');
    }
    return msaNew;
}

std::tuple<std::string, std::string> refineOne(
    int8_t * tinySubMatAA,
    int8_t * tinySubMat3Di,
    std::string msaAa,
    std::string msa3Di,
    PSSMCalculator &calculator_aa,
    MsaFilter &filter_aa,
    SubstitutionMatrix &subMat_aa,
    PSSMCalculator &calculator_3di,
    MsaFilter &filter_3di,
    SubstitutionMatrix &subMat_3di,
    StructureSmithWaterman &structureSmithWaterman,
    bool compBiasCorrection,
    bool wg,
    float filterMaxSeqId,
    float matchRatio,
    float qsc,
    float scoreBias3Di,
    float scoreBiasAa,
    int Ndiff,
    int filterMinEnable,
    int filterMsa,
    int gapExtend,
    int gapOpen,
    int maxSeqLen,
    int sequenceCnt,
    std::string qid
) {
    std::vector<bool> group(sequenceCnt);
    std::string subMSA1_aa;
    std::string subMSA2_aa;
    std::string subMSA1_3di;
    std::string subMSA2_3di;
    
    Sequence* subSequence1_aa;
    Sequence* subSequence2_aa;
    Sequence* subSequence1_3di;
    Sequence* subSequence2_3di;
    bool *maskBool1;
    bool *maskBool2;

    // create two groups using headers
    for (int j = 0; j < sequenceCnt; j++)
        group[j] = (std::rand() % 2 == 0);

    makeSubMSA(msaAa,  subMSA1_aa,  subMSA2_aa,  group); 
    makeSubMSA(msa3Di, subMSA1_3di, subMSA2_3di, group); 

    // create profile subMSA 1 -> Sequence
    //                subMSA 2 -> Sequence
    std::string subProfile1_aa = fastamsa2profile(
        subMSA1_aa, calculator_aa, filter_aa, subMat_aa, maxSeqLen,
        sequenceCnt + 1, matchRatio, filterMsa,
        compBiasCorrection,
        qid, filterMaxSeqId, Ndiff, 0,
        qsc, filterMinEnable, wg, NULL, scoreBiasAa
    );
    std::string subProfile2_aa = fastamsa2profile(
        subMSA2_aa, calculator_aa, filter_aa, subMat_aa, maxSeqLen,
        sequenceCnt + 1, matchRatio, filterMsa,
        compBiasCorrection,
        qid, filterMaxSeqId, Ndiff, 0,
        qsc, filterMinEnable, wg, NULL, scoreBiasAa
    );
    std::string mask1;
    std::string mask2;
    for (size_t k = subProfile1_aa.length() - 1; subProfile1_aa[k] != '\n'; k--)
        mask1.push_back(subProfile1_aa[k]);
    for (size_t k = subProfile2_aa.length() - 1; subProfile2_aa[k] != '\n'; k--)
        mask2.push_back(subProfile2_aa[k]);
    std::reverse(mask1.begin(), mask1.end());
    std::reverse(mask2.begin(), mask2.end());

    subProfile1_aa.erase(subProfile1_aa.length() - mask1.length() - 1);
    subProfile2_aa.erase(subProfile2_aa.length() - mask2.length() - 1);

    maskBool1 = new bool[mask1.length()];
    maskBool2 = new bool[mask2.length()];
    for (size_t k = 0; k < mask1.length(); ++k)
        maskBool1[k] = (mask1[k] == '1') ? true : false;
    for (size_t k = 0; k < mask2.length(); ++k)
        maskBool2[k] = (mask2[k] == '1') ? true : false;
    
    std::string subProfile1_3di = fastamsa2profile(
        subMSA1_3di, calculator_3di, filter_3di, subMat_3di, maxSeqLen,
        sequenceCnt + 1, matchRatio, filterMsa,
        compBiasCorrection,
        qid, filterMaxSeqId, Ndiff, 0,
        qsc, filterMinEnable, wg, maskBool1, scoreBias3Di
    );
    std::string subProfile2_3di = fastamsa2profile(
        subMSA2_3di, calculator_3di, filter_3di, subMat_3di, maxSeqLen,
        sequenceCnt + 1, matchRatio, filterMsa,
        compBiasCorrection,
        qid, filterMaxSeqId, Ndiff, 0,
        qsc, filterMinEnable, wg, maskBool2, scoreBias3Di
    );
    assert(subProfile1_aa.length() == subProfile1_3di.length());
    assert(subProfile2_aa.length() == subProfile2_3di.length());

    subSequence1_aa  = new Sequence(maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_aa,  0, false, compBiasCorrection);
    subSequence2_aa  = new Sequence(maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_aa,  0, false, compBiasCorrection);
    subSequence1_3di = new Sequence(maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_3di, 0, false, compBiasCorrection);
    subSequence2_3di = new Sequence(maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_3di, 0, false, compBiasCorrection);

    subSequence1_aa->mapSequence(0, 0, subProfile1_aa.c_str(), subProfile1_aa.length() / Sequence::PROFILE_READIN_SIZE);
    subSequence2_aa->mapSequence(1, 1, subProfile2_aa.c_str(), subProfile2_aa.length() / Sequence::PROFILE_READIN_SIZE);
    subSequence1_3di->mapSequence(0, 0, subProfile1_3di.c_str(), subProfile1_3di.length() / Sequence::PROFILE_READIN_SIZE);
    subSequence2_3di->mapSequence(1, 1, subProfile2_3di.c_str(), subProfile2_3di.length() / Sequence::PROFILE_READIN_SIZE);

    // TODO move to function in structuremsa
    float q_neff_sum = 0.0;
    float t_neff_sum = 0.0;
    for (int i = 0; i < subSequence1_aa->L; i++)
        q_neff_sum += subSequence1_aa->neffM[i];
    for (int i = 0; i < subSequence2_aa->L; i++)
        t_neff_sum += subSequence2_aa->neffM[i];
    if (q_neff_sum <= t_neff_sum) {
        std::swap(mask1, mask2);
        std::swap(subMSA1_3di, subMSA2_3di);
        std::swap(subMSA1_aa, subMSA2_aa);
        std::swap(subSequence1_3di, subSequence2_3di);
        std::swap(subSequence1_aa, subSequence2_aa);
    }

    // pairwise align
    structureSmithWaterman.ssw_init(subSequence1_aa, subSequence1_3di, tinySubMatAA, tinySubMat3Di, &subMat_aa);
    std::vector<int> map1 = maskToMapping(mask1);
    std::vector<int> map2 = maskToMapping(mask2);
    std::vector<std::vector<std::vector<int> > > neighbours; // TODO remove this?
    Matcher::result_t result = pairwiseAlignment(
        structureSmithWaterman,
        subSequence1_aa->L,
        subSequence1_aa,
        subSequence1_3di,
        subSequence2_aa,
        subSequence2_3di,
        gapOpen,
        gapExtend,
        &subMat_aa,
        &subMat_3di,
        neighbours,
        map1,
        map2
    );
    std::vector<Instruction> qBt;
    std::vector<Instruction> tBt;
    getMergeInstructions(result, map1, map2, qBt, tBt);
    std::string merged_aa  = mergeTwoMsa(subMSA1_aa,  subMSA2_aa,  result, map1, map2, qBt, tBt);
    std::string merged_3di = mergeTwoMsa(subMSA1_3di, subMSA2_3di, result, map1, map2, qBt, tBt);

    assert(merged_aa.length() == merged_3di.length());
    
    delete subSequence1_aa;
    delete subSequence2_aa;
    delete subSequence1_3di;
    delete subSequence2_3di;
    delete maskBool1;
    delete maskBool2;
    
    return std::make_tuple(merged_aa, merged_3di);
}

std::tuple<std::string, std::string> refineMany(
    int8_t * tinySubMatAA,
    int8_t * tinySubMat3Di,
    DBReader<unsigned int> &seqDbrAA,
    DBReader<unsigned int> &seqDbr3Di,
    DBReader<unsigned int> &seqDbrCA,
    std::string msaAa,
    std::string msa3Di,
    PSSMCalculator &calculator_aa,
    MsaFilter &filter_aa,
    SubstitutionMatrix &subMat_aa,
    PSSMCalculator &calculator_3di,
    MsaFilter &filter_3di,
    SubstitutionMatrix &subMat_3di,
    StructureSmithWaterman & structureSmithWaterman,
    int iterations,
    bool compBiasCorrection,
    bool wg,
    float filterMaxSeqId,
    float matchRatio,
    float qsc,
    float scoreBias3Di,
    float scoreBiasAa,
    int Ndiff,
    int filterMinEnable,
    int filterMsa,
    int gapExtend,
    int gapOpen,
    int maxSeqLen,
    int sequenceCnt,
    std::string qid,
    float pairThreshold
) {
    std::cout << "Running " << iterations << " refinement iterations\n";

    double prevLDDT = getLDDTScore(seqDbrAA, seqDbr3Di, seqDbrCA, msaAa, pairThreshold);
    std::cout << "Initial LDDT: " << prevLDDT << '\n';

    for (int i = 0; i < iterations; i++) {
        std::tuple<std::string, std::string> msas = refineOne(
            tinySubMatAA,
            tinySubMat3Di,
            msaAa,
            msa3Di,
            calculator_aa,
            filter_aa,
            subMat_aa,
            calculator_3di,
            filter_3di,
            subMat_3di,
            structureSmithWaterman,
            compBiasCorrection,
            wg,
            filterMaxSeqId,
            matchRatio,
            qsc,
            scoreBias3Di,
            scoreBiasAa,
            Ndiff,
            filterMinEnable,
            filterMsa,
            gapExtend,
            gapOpen,
            maxSeqLen,
            sequenceCnt,
            qid
        );
        float lddtScore = getLDDTScore(seqDbrAA, seqDbr3Di, seqDbrCA, std::get<0>(msas), pairThreshold);
        if (lddtScore > prevLDDT) {
            std::cout << std::fixed << std::setprecision(4) << prevLDDT << " -> " << lddtScore << " (+" << (lddtScore - prevLDDT) << ") #" << i << '\n';
            prevLDDT = lddtScore;
            msaAa  = std::get<0>(msas);
            msa3Di = std::get<1>(msas);
        }
    }
    return std::make_tuple(msaAa, msa3Di);
}

std::pair<std::string, int> parseMSADb(std::string db) {
    std::string msa;
    int count = 0;
    KSeqWrapper* kseq = KSeqFactory(db.c_str());
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;
        msa.push_back('>');
        msa.append(entry.name.s); 
        msa.push_back('\n');
        msa.append(entry.sequence.s);
        msa.push_back('\n');
        count++;
    }
    return std::make_pair(msa, count);
}

int refinemsa(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    
    if (par.refineIters <= 0.0) {
        std::cerr << "Expected >0 refinement iterations\n";
        return EXIT_FAILURE;
    }

    DBReader<unsigned int> seqDbrAA(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_LOOKUP_REV);
    seqDbrAA.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbr3Di((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbr3Di.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbrCA((par.db1+"_ca").c_str(), (par.db1+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbrCA.open(DBReader<unsigned int>::NOSORT);

    IndexReader qdbrH(par.db1, par.threads, IndexReader::HEADERS, touch ? IndexReader::PRELOAD_INDEX : 0);

    // Read in FASTA alignment
    // Also get # structures here for sequenceCnt
    // Use AA MSA to generate 3Di ones
    int sequenceCnt;
    std::string msaAa;
    std::tie(msaAa, sequenceCnt) = parseMSADb(par.db2);
    std::string msa3Di = convertMSA(seqDbrAA, seqDbr3Di, msaAa);
    
    SubstitutionMatrix subMat_3di(par.scoringMatrixFile.values.aminoacid().c_str(), par.bitFactor3Di, par.scoreBias3di);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
        }
    }
    SubstitutionMatrix subMat_aa(blosum.c_str(), par.bitFactorAa, par.scoreBiasAa);

    // Substitution matrices needed for query profile
    int8_t *tinySubMatAA  = (int8_t*) mem_align(ALIGN_INT, subMat_aa.alphabetSize * 32);
    int8_t *tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat_3di.alphabetSize * 32);

    for (int i = 0; i < subMat_aa.alphabetSize; i++)
        for (int j = 0; j < subMat_aa.alphabetSize; j++)
            tinySubMatAA[i * subMat_aa.alphabetSize + j] = subMat_aa.subMatrix[i][j];
    for (int i = 0; i < subMat_3di.alphabetSize; i++)
        for (int j = 0; j < subMat_3di.alphabetSize; j++)
            tinySubMat3Di[i * subMat_3di.alphabetSize + j] = subMat_3di.subMatrix[i][j]; // for farrar profile

    StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat_3di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMat_aa, &subMat_3di);
    MsaFilter filter_aa(par.maxSeqLen + 1, sequenceCnt + 1, &subMat_aa, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
    MsaFilter filter_3di(par.maxSeqLen + 1, sequenceCnt + 1, &subMat_3di, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid()); 
    PSSMCalculator calculator_aa(&subMat_aa, par.maxSeqLen + 1, sequenceCnt + 1, par.pcmode, par.pcaAa, par.pcbAa
#ifdef GAP_POS_SCORING
            , par.gapOpen.values.aminoacid(), par.gapPseudoCount
#endif
            );
    PSSMCalculator calculator_3di(&subMat_3di, par.maxSeqLen + 1, sequenceCnt + 1, par.pcmode, par.pca3di, par.pcb3di
#ifdef GAP_POS_SCORING
            , par.gapOpen.values.aminoacid(), par.gapPseudoCount
#endif
            );
    
    // Refine for N iterations
    // std::tuple<std::string, std::string, std::string> newMSA = refineMany(
    std::tie(msaAa, msa3Di) = refineMany(
        tinySubMatAA,
        tinySubMat3Di,
        seqDbrAA,
        seqDbr3Di,
        seqDbrCA,
        msaAa,
        msa3Di,
        calculator_aa,
        filter_aa,
        subMat_aa,
        calculator_3di,
        filter_3di,
        subMat_3di,
        structureSmithWaterman,
        par.refineIters,
        par.compBiasCorrection,
        par.wg,
        par.filterMaxSeqId,
        par.matchRatio,
        par.qsc,
        par.scoreBias3di,
        par.scoreBiasAa,
        par.Ndiff,
        par.filterMinEnable,
        par.filterMsa,
        par.gapExtend.values.aminoacid(),
        par.gapOpen.values.aminoacid(),
        par.maxSeqLen,
        sequenceCnt,
        par.qid,
        par.pairThreshold
    );
    
    // Write final MSA to file
    DBWriter resultWriter(
        par.db3.c_str(),
        (par.db3 + ".index").c_str(),
        static_cast<unsigned int>(par.threads),
        par.compressed,
        Parameters::DBTYPE_OMIT_FILE
    );
    resultWriter.open();
    std::string buffer;
    KSeqBuffer *kseq = new KSeqBuffer(msaAa.c_str(), msaAa.length());
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;
        buffer.append(1, '>');
        buffer.append(entry.name.s);
        buffer.append(1, '\n');
        buffer.append(entry.sequence.s, entry.sequence.l);
        buffer.append(1, '\n');
        resultWriter.writeAdd(buffer.c_str(), buffer.size(), 0);
        buffer.clear();
    }
    resultWriter.writeEnd(0, 0, false, 0);
    resultWriter.close(true);
    FileUtil::remove((par.db3 + ".index").c_str());

    // Cleanup
    seqDbrAA.close();
    seqDbr3Di.close();
    delete [] tinySubMatAA;
    delete [] tinySubMat3Di;

    return EXIT_SUCCESS;
}