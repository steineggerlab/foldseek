#include "DBReader.h"
#include "IndexReader.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Parameters.h"
#include "StructureUtil.h"
#include "Util.h"
#include <cassert>

#include "kseq.h"
#include "KSeqBufferReader.h"
#include "KSeqWrapper.h"
#include "LDDT.h"

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
    int sequenceCnt = seqDbr3Di.getSize();
    std::vector<std::string> headers(sequenceCnt);
    std::vector<std::string> sequences(sequenceCnt);
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;
        size_t idx = seqDbrAA.getLookupIdByAccession(entry.name.s);
        std::cout << entry.name.s << " " << idx << std::endl;
        headers[idx] = entry.name.s;
        sequences[idx] = entry.sequence.s;
        if (alnLength == 0)
            alnLength = (int)entry.sequence.l;
    }

    std::vector<bool> ids(sequenceCnt);  
    std::fill(ids.end() - 2, ids.end(), true);
    double sum = 0.0;
    int numPairs = 0;
    do {
        numPairs++;
        
        // Get the qId/tId of this permutation
        int qId = -1;
        int tId = -1;
        for (int i = 0; i < sequenceCnt; ++i) {
            if (ids[i] && qId == -1) { qId = i; continue; }
            if (ids[i] && tId == -1) { tId = i; break; }
        }

        // Generate a CIGAR string from qId-tId sub-alignment, ignoring -- columns
        // e.g. --X-XX-X---XX-
        //      Y---YYYY---YYY
        //          MMDM   MM
        Matcher::result_t res;
        res.backtrace = "";
        bool started = false;  // flag for first Match column in backtrace
        int q  = 0;  // index in msa
        int qr = 0;  // index in seq
        int t  = 0;
        int tr = 0;
        while (q < alnLength && t < alnLength) {
            char qc = sequences[qId][q];
            char tc = sequences[tId][t];
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
                res.backtrace.push_back('M');
                ++qr;
                ++tr;
                res.qEndPos = qr;
                res.dbEndPos = tr;
            }
            ++q;
            ++t;
        }

        // These should be equal to qr/tr 
        int qLen = seqDbrAA.getSeqLen(qId);
        int tLen = seqDbr3Di.getSeqLen(tId);
        
        // Remove D/I from backtrace after last M
        res.backtrace.shrink_to_fit();
        size_t i;
        for (i = res.backtrace.length() - 1; res.backtrace[i] != 'M'; i--);
        res.backtrace.erase(i + 1, res.backtrace.length() - i);

        // Get c-alpha info from _ca db for q and t
        float *queryCaData  = (float*)seqDbrCA.sequenceReader->getData(qId, 0);
        float *targetCaData = (float*)seqDbrCA.sequenceReader->getData(tId, 0);
        
        // Calculate LDDT using created result_t object
        LDDTCalculator *lddtcalculator = new LDDTCalculator(qLen, tLen);
        lddtcalculator->initQuery(qLen, queryCaData, &queryCaData[qLen], &queryCaData[2 * qLen]);
        LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(
            tLen,
            res.qStartPos,
            res.dbStartPos,
            res.backtrace,
            targetCaData,
            &targetCaData[tLen],
            &targetCaData[2 * tLen]
        );
        
        // TODO: output per-column LDDT
        // e.g. X   X   -   X
        //      Y   Y   Y   Y
        //      Z   Z   -   Z
        //     0.4 0.7  0  0.6
        // - Map perCaLddtScore values to MSA column indices
        
        std::cout << headers[qId] << " " << headers[tId] << " " << lddtres.avgLddtScore << std::endl;
        // if (!std::isnan(lddtres.avgLddtScore))
        sum += lddtres.avgLddtScore;

        delete lddtcalculator;
    } while (std::next_permutation(ids.begin(), ids.end()));
    
    std::cout << "Average MSA LDDT: " << sum / (double)numPairs << std::endl;
    return EXIT_SUCCESS;
}