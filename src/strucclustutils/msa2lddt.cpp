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

// LDDT scoring
// - Iterate unique pairs of sequences in MSA
//   - Generate faux Matcher::result_t
//   - Initialise LDDTCalculator using q/t seq len
//   - Calculate LDDT
//   - Add to sum
int msa2lddt(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    DBReader<unsigned int> seqDbrAA(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_LOOKUP_REV);
    seqDbrAA.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbr3Di((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbr3Di.open(DBReader<unsigned int>::NOSORT);
    IndexReader seqDbrCA(par.db1, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY), touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca");
    IndexReader headerDB(par.db1, par.threads, IndexReader::HEADERS, touch ? IndexReader::PRELOAD_INDEX : 0);
    
    // Read in MSA, mapping headers to database indices
    KSeqWrapper* kseq = KSeqFactory(par.db2.c_str());
    int alnLength = 0;
    int numSequences = 0;
    int maxSeqLength = 0;
    int sequenceCnt = seqDbr3Di.getSize();
    std::vector<std::string> headers;
    std::vector<size_t> indices;
    std::vector<std::string> sequences;
    std::vector<std::string> sequences3di;
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;
        size_t idx = seqDbrAA.getLookupIdByAccession(entry.name.s);
        unsigned int key = seqDbrAA.getLookupKey(idx);
        size_t seqIdAA = seqDbrAA.getId(key);
        if (seqIdAA == UINT_MAX)
            Debug(Debug::WARNING) << "Key not found in seqDbrAA: " << key << "\n";
        size_t seqId3Di = seqDbr3Di.getId(key);
        if (seqId3Di == UINT_MAX)
            Debug(Debug::WARNING) << "Key not found in seqDbr3di: " << key << "\n";
        headers.push_back(entry.name.s);
        indices.push_back(seqIdAA);
        std::string seqAA = seqDbrAA.getData(seqIdAA, 0);
        seqAA.pop_back();
        std::string seq3Di = seqDbr3Di.getData(seqId3Di, 0);
        seq3Di.pop_back();
        for (size_t i = 0; i < entry.sequence.l; i++) {
            if (entry.sequence.s[i] == '-') {
                seqAA.insert(i, "-");
                seq3Di.insert(i, "-");
            }
        }
        sequences.push_back(seqAA);
        sequences.push_back(seq3Di);
        if (alnLength == 0)
            alnLength = (int)entry.sequence.l;
        if (seqAA.length() > maxSeqLength)
            maxSeqLength = seqAA.length();
        numSequences++;
    }
    
    // Track per-column scores and no. non-gaps to avg
    std::vector<float> perColumnScore(alnLength, 0.0);
    std::vector<float> perColumnSOP(alnLength, 0.0);
    std::vector<float> perColumnCount(alnLength, 0.0);

    std::vector<bool> ids(numSequences);  
    std::fill(ids.end() - 2, ids.end(), true);
    double sum = 0.0;
    int numPairs = 0;

    LDDTCalculator *lddtcalculator = new LDDTCalculator(maxSeqLength + 1, maxSeqLength + 1);
    int prevqId = -1;
    do {
        numPairs++;
        
        // Get the qId/tId of this permutation
        int qId = -1;
        int tId = -1;
        for (int i = 0; i < sequenceCnt; ++i) {
            if (ids[i] && qId == -1) { qId = i; continue; }
            if (ids[i] && tId == -1) { tId = i; break; }
        }
        
        // lengths
        // per column lddt idx <-> msa idx
        std::vector<int> match_to_msa;
        
        // Generate a CIGAR string from qId-tId sub-alignment, ignoring -- columns
        // e.g. --X-XX-X---XX-
        //      Y---YYYY---YYY
        //          MMDM   MM
        Matcher::result_t res;
        res.backtrace = "";
        bool started = false;  // flag for first Match column in backtrace
        int m  = 0;  // index in msa
        int qr = 0;  // index in seq
        int tr = 0;
        while (m < alnLength) {
            char qc = sequences[qId][m];
            char tc = sequences[tId][m];
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
        size_t i;
        for (i = res.backtrace.length() - 1; res.backtrace[i] != 'M'; i--);
        res.backtrace.erase(i + 1, res.backtrace.length() - i);

        // Get c-alpha info from _ca db for q and t
        Coordinate16 qcoords;
        char *qcadata = seqDbrCA.sequenceReader->getData(indices[qId], 0);
        size_t qCaLength = seqDbrCA.sequenceReader->getEntryLen(indices[qId]);
        float *queryCaData = qcoords.read(qcadata, qr, qCaLength);

        Coordinate16 tcoords;
        char *tcadata = seqDbrCA.sequenceReader->getData(indices[tId], 0);
        size_t tCaLength = seqDbrCA.sequenceReader->getEntryLen(indices[tId]);
        float *targetCaData = tcoords.read(tcadata, tr, tCaLength);

        // Calculate LDDT using created result_t object
        // Only re-initialise when query actually changes
        if (prevqId != qId)
            lddtcalculator->initQuery(qr, queryCaData, &queryCaData[qr], &queryCaData[2 * qr]);
        prevqId = qId;

        LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(
            tr,
            res.qStartPos,
            res.dbStartPos,
            res.backtrace,
            targetCaData,
            &targetCaData[tr],
            &targetCaData[2 * tr]
        );
        
        for (int i = 0; i < lddtres.scoreLength; i++) {
            if (lddtres.perCaLddtScore[i] == 0.0) continue;
            size_t idx = match_to_msa[i];
            perColumnSOP[idx] += (sequences[qId][idx] == sequences[tId][idx]);
            perColumnScore[idx] += lddtres.perCaLddtScore[i];
            perColumnCount[idx]++;
        }
        
        sum += lddtres.avgLddtScore;
    } while (std::next_permutation(ids.begin(), ids.end()));
    delete lddtcalculator;
   
    std::cout << "Average MSA LDDT: " << sum / (double)numPairs << std::endl;
    
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
            if (perColumnCount[i] == 0) continue;
            std::string entry;
            entry.append("[");
            entry.append(std::to_string(i));
            entry.append(",");
            entry.append(std::to_string(perColumnScore[i] / (double)perColumnCount[i]));
            entry.append(",");
            entry.append(std::to_string(perColumnSOP[i] / (double)perColumnCount[i]));
            entry.append("],");
            resultWriter.writeData(entry.c_str(), entry.length(), 0, 0, false, false);
        }
        
        std::string end = "],";
        end.append("{'db':'");
        end.append(par.db1);
        end.append("','msa':'");
        end.append(par.db2);
        end.append("','lddt':");
        end.append(std::to_string(sum / numPairs));
        end.append("});</script>");

        resultWriter.writeData(end.c_str(), end.length(), 0, 0, false, false);
        resultWriter.writeEnd(0, 0, false, 0);
        resultWriter.close(true);
        FileUtil::remove(lddtHtmlIdx.c_str());
    }

    return EXIT_SUCCESS;
}