#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Parameters.h"
#include "StructureUtil.h"
#include "Util.h"
#include <cassert>

#define ZSTD_STATIC_LINKING_ONLY

#include <zstd.h>
#include "msa.html.zst.h"

#include "kseq.h"
#include "KSeqBufferReader.h"
#include "KSeqWrapper.h"
#include "LDDT.h"
#include "Coordinate16.h"

#ifdef OPENMP
#include <omp.h>
#endif

#pragma omp declare reduction(vsum : std::vector<float> : \
    std::transform(omp_in.begin(), omp_in.end(), \
                   omp_out.begin(), omp_out.begin(), std::plus<float>())) \
    initializer(omp_priv(omp_orig))

#pragma omp declare reduction(vsum : std::vector<int> : \
    std::transform(omp_in.begin(), omp_in.end(), \
                   omp_out.begin(), omp_out.begin(), std::plus<int>())) \
    initializer(omp_priv(omp_orig))


std::tuple<std::vector<float>, std::vector<int>, float> calculate_lddt(
    std::vector<std::string> &sequences,
    std::vector<size_t> &indices,
    std::vector<int> &lengths,  // ungapped
    DBReader<unsigned int> * seqDbrCA,
    float pairThreshold
) {
    // Track per-column scores and no. non-gaps to avg
    std::vector<float> perColumnScore(sequences[0].length(), 0.0);
    std::vector<int>   perColumnCount(sequences[0].length(), 0);

    float sum = 0.0;
    int numPairs = sequences.size() * (sequences.size() - 1) / 2;
    size_t maxSeqLength = 0;
    for (size_t i = 0; i < sequences.size(); i++) {
        if (sequences[i].length() > maxSeqLength)
            maxSeqLength = sequences[i].length();
    }
    
#pragma omp parallel reduction(+:sum) reduction(vsum:perColumnScore,perColumnCount)
{
    LDDTCalculator *lddtcalculator = new LDDTCalculator(maxSeqLength + 1, maxSeqLength + 1);
    std::vector<Matcher::result_t> results(numPairs);
    
    for (size_t i = 0; i < sequences.size(); i++) {
        Coordinate16 qcoords;
        char *qcadata = seqDbrCA->getData(indices[i], 0);
        size_t qCaLength = seqDbrCA->getEntryLen(indices[i]);
        float *queryCaData = qcoords.read(qcadata, lengths[i], qCaLength);
        
        lddtcalculator->initQuery(
            lengths[i],
            queryCaData,
            &queryCaData[lengths[i]],
            &queryCaData[lengths[i] * 2]
        );

#pragma omp for schedule(dynamic, 1)
        for (size_t j = i + 1; j < sequences.size(); j++) {
            Coordinate16 tcoords;
            char *tcadata = seqDbrCA->getData(indices[j], 0);
            size_t tCaLength = seqDbrCA->getEntryLen(indices[j]);
            float *targetCaData = tcoords.read(tcadata, lengths[j], tCaLength);

            // lengths per column lddt idx <-> msa idx
            std::vector<int> match_to_msa;

            // Generate a CIGAR string from qId-tId sub-alignment, ignoring -- columns
            // e.g. --X-XX-X---XX-
            //      Y---YYYY---YYY
            //          MMDM   MM
            Matcher::result_t res;
            res.backtrace = "";
            bool started = false;  // flag for first Match column in backtrace
            size_t m  = 0;  // index in msa
            size_t qr = 0;  // index in seq
            size_t tr = 0;
            while (m < sequences[0].length()) {
                char qc = sequences[i][m];
                char tc = sequences[j][m];
                if (qc == '-' && tc == '-') {
                    // Skip gap columns
                } else if (qc == '-') {
                    if (started) res.backtrace.push_back('D');
                    ++tr;
                } else if (tc == '-') {
                    if (started) res.backtrace.push_back('I');
                    ++qr;
                } else {
                    if (!started) {
                        started = true; 
                        res.qStartPos = qr;
                        res.dbStartPos = tr;
                    }
                    match_to_msa.push_back(m);
                    res.backtrace.push_back('M');
                    res.qEndPos = qr;
                    res.dbEndPos = tr;
                    ++qr;
                    ++tr;
                }
                ++m;
            }
            
            // If no alignment between the two sequences, skip
            if (res.backtrace.length() == 0)
                continue;

            // Remove D/I from backtrace after last M
            res.backtrace.shrink_to_fit();
            size_t k;
            for (k = res.backtrace.length() - 1; res.backtrace[k] != 'M'; k--);
            res.backtrace.erase(k + 1);

            LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(
                lengths[j],
                res.qStartPos,
                res.dbStartPos,
                res.backtrace,
                targetCaData,
                &targetCaData[lengths[j]],
                &targetCaData[lengths[j] * 2]
            );
            
            if (std::isnan(lddtres.avgLddtScore)) {
                std::cout << "Found nan\n";
                // std::cout << sequences[i]  << '\n';
                // std::cout << res.backtrace << '\n';
                // std::cout << sequences[j]  << '\n';
                // for (int k = 0; k < lddtres.scoreLength; k++) {
                //     std::cout << lddtres.perCaLddtScore[k] << " ";
                // }
                // std::cout << '\n';
            }
            
            for (int k = 0; k < lddtres.scoreLength; k++) {
                if (lddtres.perCaLddtScore[k] == 0.0)
                    continue;
                int idx = match_to_msa[k];
                perColumnCount[idx] += 1;
                perColumnScore[idx] += lddtres.perCaLddtScore[k];
            }
            sum += lddtres.avgLddtScore;
        }
    }
    delete lddtcalculator;
}

    float scaledSum = 0.0;
    int numCols = 0;
    for (size_t i = 0; i < perColumnCount.size(); i++) {
        float pairSupport = perColumnCount[i] / static_cast<float>(numPairs);

        float residuesInColumn = 0.0;
        for (size_t j = 0; j < sequences.size(); j++) {
            if (sequences[j][i] != '-') {
                residuesInColumn++;
            }
        }
        // TODO maybe change to bool flag to just count/not count single-residue columns in entire MSA score
        if ((residuesInColumn / sequences.size()) < pairThreshold) {
            perColumnScore[i] = 0.0;
            perColumnCount[i] = 0;
            continue;
        }
        if (perColumnCount[i] > 0) {
            perColumnScore[i] /= perColumnCount[i];                                   // get mean LDDT for this column
            // perColumnScore[i] *= pairSupport;  // scale by % of pairs with LDDT >0
        } else {
            perColumnScore[i] = 0.0;
        }
        // scaledSum += perColumnScore[i];
        // if (perColumnCount[i] / static_cast<float>(numPairs) >= pairThreshold) {
        scaledSum += perColumnScore[i];
        numCols++;
        // }
    }
    // float lddtScore = sum / static_cast<double>(numPairs);
    // float lddtScore = (scaledSum / perColumnCount.size());  // get mean over all columns
    float lddtScore = scaledSum / static_cast<float>(numCols);                                        

    return std::make_tuple(perColumnScore, perColumnCount, lddtScore);
}

void parseFasta(
    KSeqWrapper *kseq,
    DBReader<unsigned int> * seqDbrAA,
    DBReader<unsigned int> * seqDbr3Di,
    std::vector<std::string> &headers,
    std::vector<size_t>      &indices,
    std::vector<int>         &lengths,
    std::vector<std::string> &sequences,
    std::vector<std::string> &sequences3di,
    int &alnLength
) {
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;
        size_t idx = seqDbrAA->getLookupIdByAccession(entry.name.s);
        unsigned int key = seqDbrAA->getLookupKey(idx);
        size_t seqIdAA = seqDbrAA->getId(key);
        if (seqIdAA == UINT_MAX)
            Debug(Debug::WARNING) << "Key not found in seqDbrAA: " << key << "\n";
        size_t seqId3Di = seqDbr3Di->getId(key);
        if (seqId3Di == UINT_MAX)
            Debug(Debug::WARNING) << "Key not found in seqDbr3di: " << key << "\n";
        headers.push_back(entry.name.s);
        indices.push_back(seqIdAA);
        std::string seqAA = seqDbrAA->getData(seqIdAA, 0);
        seqAA.pop_back();
        std::string seq3Di = seqDbr3Di->getData(seqId3Di, 0);
        seq3Di.pop_back();
        lengths.push_back(seqAA.length());
        for (size_t i = 0; i < entry.sequence.l; i++) {
            if (entry.sequence.s[i] == '-') {
                seqAA.insert(i, "-");
                seq3Di.insert(i, "-");
            }
        }
        sequences.push_back(seqAA);
        sequences3di.push_back(seq3Di);
        if (alnLength == 0)
            alnLength = (int)entry.sequence.l;
    }
}

float getLDDTScore(
    DBReader<unsigned int> &seqDbrAA,
    DBReader<unsigned int> &seqDbr3Di,
    DBReader<unsigned int> &seqDbrCA,
    std::string msa,
    float pairThreshold
) {
    KSeqWrapper* kseq = new KSeqBuffer(msa.c_str(), msa.length());
    int alnLength = 0;
    std::vector<std::string> hdrs;
    std::vector<size_t>      inds;
    std::vector<int>         lens;
    std::vector<std::string> seqs;
    std::vector<std::string> seqs3di;
    parseFasta(kseq, &seqDbrAA, &seqDbr3Di, hdrs, inds, lens, seqs, seqs3di, alnLength);
    delete kseq;

    std::vector<float> perColumnScore;
    std::vector<int>   perColumnCount;
    float lddtScore;
    std::tie(perColumnScore, perColumnCount, lddtScore) = calculate_lddt(seqs, inds, lens, &seqDbrCA, pairThreshold);

    return lddtScore;
}

int msa2lddt(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    DBReader<unsigned int> seqDbrAA(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_LOOKUP_REV);
    seqDbrAA.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbr3Di((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbr3Di.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbrCA((par.db1+"_ca").c_str(), (par.db1+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbrCA.open(DBReader<unsigned int>::NOSORT);
    IndexReader headerDB(par.db1, par.threads, IndexReader::HEADERS, touch ? IndexReader::PRELOAD_INDEX : 0);
    
    // Read in MSA, mapping headers to database indices
    KSeqWrapper* kseq = KSeqFactory(par.db2.c_str());
    int alnLength = 0;
    std::vector<std::string> headers;
    std::vector<size_t> indices;
    std::vector<int> lengths;
    std::vector<std::string> sequences;
    std::vector<std::string> sequences3di;
    parseFasta(kseq, &seqDbrAA, &seqDbr3Di, headers, indices, lengths, sequences, sequences3di, alnLength);

    // Calculate LDDT
    std::vector<float> perColumnScore;
    std::vector<int>   perColumnCount;
    float lddtScore;
    std::tie(perColumnScore, perColumnCount, lddtScore) = calculate_lddt(sequences, indices, lengths, &seqDbrCA, par.pairThreshold);
    
    std::cout << "Average MSA LDDT: " << lddtScore << std::endl;
    
    // Write clustal format MSA HTML
    // TODO: make optional
    if (par.lddtHtml != "") {
        std::string lddtHtmlIdx = par.lddtHtml + ".index";
        DBWriter resultWriter(par.lddtHtml.c_str(), lddtHtmlIdx.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_OMIT_FILE);
        resultWriter.open();
        
        // Read in template and write to .html
        size_t dstSize = ZSTD_findDecompressedSize(msa_html_zst, msa_html_zst_len);
        char* dst = (char*)malloc(sizeof(char) * dstSize);
        size_t realSize = ZSTD_decompress(dst, dstSize, msa_html_zst, msa_html_zst_len);
        resultWriter.writeData(dst, realSize, 0, 0, false, false);

        // Aligned sequences, as [[header, sequence], [header, sequence], ...]
        const char* scriptBlock = "<script>render([";
        resultWriter.writeData(scriptBlock, strlen(scriptBlock), 0, 0, false, false);
        for (size_t i = 0; i < sequences.size(); i++) {
            std::string entry;
            entry.append("['");
            entry.append(headers[i]);
            entry.append("','");
            entry.append(sequences[i]);
            entry.append("','");
            entry.append(sequences3di[i]);
            entry.append("'],");
            resultWriter.writeData(entry.c_str(), entry.length(), 0, 0, false, false);
        }
        std::string middle = "],[";
        resultWriter.writeData(middle.c_str(), middle.length(), 0, 0, false, false);

        // Per-column scores, as [[id, score], [id, score], ...]
        // TODO: optionally save this as .csv
        for (int i = 0; i < alnLength; i++) {
            if (perColumnCount[i] == 0)
                continue;
            std::string entry;
            entry.append("[");
            entry.append(std::to_string(i));
            entry.append(",");
            entry.append(std::to_string(perColumnScore[i]));   // / (double)perColumnCount[i]));
            entry.append("],");
            resultWriter.writeData(entry.c_str(), entry.length(), 0, 0, false, false);
        }
        
        std::string end = "],";
        end.append("{'db':'");
        end.append(par.db1);
        end.append("','msa':'");
        end.append(par.db2);
        end.append("','lddt':");
        end.append(std::to_string(lddtScore));
        end.append("});</script>");

        resultWriter.writeData(end.c_str(), end.length(), 0, 0, false, false);
        resultWriter.writeEnd(0, 0, false, 0);
        resultWriter.close(true);
        FileUtil::remove(lddtHtmlIdx.c_str());
    }

    return EXIT_SUCCESS;
}