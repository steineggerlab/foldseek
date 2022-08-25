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
#include "BacktraceTranslator.h"
#include "MultipleAlignment.h"
#include "Sequence.h"
#include <tuple>
#include <set>

#include "MSANode.h"

#ifdef OPENMP
#include <omp.h>
#endif

#define	EXIT_FAILURE	1
#define	EXIT_SUCCESS	0

// Tests c is a match
inline bool match(char c) {
    return (c >= 'A' && c <= 'Z') || c == '-';
}

inline bool insertion(char c) {
    return (c >= 'a' && c <= 'z') || c == '.';
}

// Transforms chr into an uppercase character
inline char uprchr(char chr) {
  return (chr >= 'a' && chr <= 'z') ? chr + 'A' - 'a' : chr;
}

// Transforms chr into lowercase character
inline char lwrchr(char chr) {
    return (chr >= 'A' && chr <= 'Z') ? chr - 'A' + 'a' : chr;
}

struct CompareResult {
    public:
        template <class T> bool operator() (T &a, T &b) {
            return a.second.score > b.second.score;
        }
} BestResult;

template <class T> void printVector(std::vector<T> vector) {
    for (size_t i = 0; i < vector.size(); i++) {
        std::cout << vector[i] << " ";
        if (i == vector.size() - 1) {
            std::cout << std::endl;
        }
    }
}

/**
 * Print two sequences, shifted by alignment result position
 */
void print_sequences(std::string a, std::string b, Matcher::result_t result) {
    if (result.qStartPos > result.dbStartPos) {
        std::cout << a;
        std::cout << std::string(result.qStartPos - result.dbStartPos, ' ') << b << std::endl;
    } else {
        std::cout << std::string(result.dbStartPos - result.qStartPos, ' ') << a;
        std::cout << b << std::endl;
    }
}

void printCurSeq(char cur_seq[], size_t length) {
    std::cout << "cur_seq: ";
    for (size_t i = 0; i < length; i++)
        std::cout << cur_seq[i];
    std::cout << std::endl;
}

/**
 * Merge alignment B into alignment A.
 */
void merge_a3m(std::string _a, std::string _b, Matcher::result_t result) {
    std::string a = _a;
    std::string b = _b;

    size_t i, j;
    std::cout << "        Score: " << result.seqId * 100 << " " << result.eval << std::endl;
    std::cout << "  qStart/qEnd: " << result.qStartPos << " " << result.qEndPos << std::endl;
    std::cout << "  tStart/tEnd: " << result.dbStartPos << " " << result.dbEndPos << std::endl;
    std::cout << "    alnLength: " << result.alnLength << std::endl;
    std::cout << "        CIGAR: " << result.backtrace << " (" << result.backtrace.length() << ")" << std::endl;

    // CIGAR backtrace relates query to target
    // i.e. 'D' --> deletion in query
    //      'I' --> insertion in query, so gap in target
    // NOTE: alnLength != CIGAR backtrace length
    size_t aIdx = result.qStartPos;
    size_t bIdx = result.dbStartPos;
    for (i = 0; i < result.backtrace.length(); i++) {
        switch (result.backtrace[i]) {
            case 'M':
                aIdx++;
                bIdx++;
                break;
            case 'D':
                a.insert(aIdx++, 1, '-');
                bIdx++;
                break;
            case 'I':
                b.insert(bIdx++, 1, '-');
                aIdx++;
                break;
        }
    }

    // Regular MSA style
    std::string diff(result.backtrace.length(), ' ');
    for (i = 0; i < result.backtrace.length(); i++) {
        if (a[i + result.qStartPos] == b[i + result.dbStartPos])
            diff[i] =  a[i + result.qStartPos];
    }
    std::cout << "\n        Query: ";
    for (i = 0; i < result.backtrace.length(); i++) {
        std::cout << a[i + result.qStartPos];
    }
    std::cout << "\n         Diff: ";
    for (i = 0; i < result.backtrace.length(); i++) {
        std::cout << diff[i];
    }
    std::cout << "\n       Target: ";
    for (i = 0; i < result.backtrace.length(); i++) {
        std::cout << b[i + result.dbStartPos];
    }

    // A3M style
    std::string fasta_a = a;
    std::string fasta_b = b;
    a = _a;
    b = _b;
    for (aIdx = result.qStartPos, bIdx = result.dbStartPos, i = 0; i < result.backtrace.length(); i++) {
        switch (result.backtrace[i]) {
            case 'M':
                aIdx++;
                bIdx++;
                break;
            case 'D':
                b[bIdx] = lwrchr(b[bIdx]);
                bIdx++;
                a.insert(aIdx++, 1, '.');
                break;
            case 'I':
                b.insert(bIdx++, 1, '-');
                aIdx++;
                break;
        }
    }
    std::cout << "\n\n        Query: ";
    for (i = 0; i < result.backtrace.length(); i++) {
        std::cout << a[i + result.qStartPos];
    }
    std::cout << "\n       Target: ";
    for (i = 0; i < result.backtrace.length(); i++) {
        std::cout << b[i + result.dbStartPos];
    }
    a.insert(0, result.dbStartPos, '-');
    b.insert(0, result.qStartPos, '-');

    std::string a3m_a = a;
    std::string a3m_b = b;

    /* int *imatch = new int[result.dbEndPos + 1]; */
    a = _a;
    b = _b;
    std::map<int, int> imap;

    // Use CIGAR string to map matching target & query indices
    std::cout << std::endl;
    for (aIdx = result.qStartPos, bIdx = result.dbStartPos, i = 0; i < result.backtrace.length(); i++) {
        switch (result.backtrace[i]) {
            case 'M':
                aIdx++;
                bIdx++;
                break;
            case 'D':
                bIdx++;
                break;
            case 'I':
                aIdx++;
                break;
        }
        imap[bIdx] = aIdx;
    }

    int par_maxcol = 32000;
    char *cur_seq = new char[par_maxcol];

    char c;
    size_t h;
    size_t l, ll;

    // Add left-end gaps
    for (h = 0; h < (size_t)result.qStartPos; h++)
        cur_seq[h] = '-';

    // Mark insertions until first match position
    j = 0;
    if (result.dbStartPos > 0) {
        for (; j < (size_t)result.dbStartPos; j++) {
            c = b[j];
            cur_seq[h++] = (c == '-') ? '.' : lwrchr(c);
        }
    }

    // If the first position is an insertion, add lowercase until match
    c = b[j];
    while (insertion(c) && c != '\0') {
        cur_seq[h++] = lwrchr(c);
        c = b[j++];
    }

    // Get index of first match state in backtrace
    l = result.backtrace.find_first_of('M');

    // Write first match state to cur_seq
    int iprev = result.qStartPos;
    int lprev = result.dbStartPos + l;

    cur_seq[h++] = b[lprev];        // First column is Match-Match state

    // Iterate next match states
    for (j = result.dbStartPos + 1; j < (result.dbStartPos + result.backtrace.length()); j++) {
        // Get index in query from imap
        i = imap[j];

        // Advance to position of next T match state j
        while ((c = b[++l]) > '\0' && insertion(c));

        int di = i - iprev;
        int dl = l - lprev;

        std::cout
            << "j: " << j << ", i: " << i << ", iprev: " << iprev << ", di: " << di
            << ", l: " << l << ", lprev: " << lprev << ", dl: " << dl
            << std::endl;
        
        if (di == 1) {
            for (ll = lprev + 1; ll < l; ll++)
                if (b[ll] != '-' && b[ll] != '.')
                    cur_seq[h++] = lwrchr(b[ll]);
            cur_seq[h++] = b[ll++];

        } else if (di == 0) {
            for (ll = lprev + 1; ll <= l; ll++)
                if (b[ll] != '-' && b[ll] != '.')
                    cur_seq[h++] = lwrchr(b[ll]);

        } else if (di >= dl) {
            for (ll = lprev + 1; ll <= (size_t)(lprev + (int)(dl / 2)); ll++)
                cur_seq[h++] = uprchr(b[ll]);
            for (int gap = 1; gap <= di - dl; gap++)
                cur_seq[h++] = '-';
            for (; ll <= l; ll++)
                cur_seq[h++] = uprchr(b[ll]);

        } else if (di < dl) {
            for (ll = lprev + 1; ll <= (size_t)(lprev + (int)(di / 2)); ll++)
                cur_seq[h++] = uprchr(b[ll]);
            for (int ins = 1; ins <= dl - di; ins++, ll++)
                if (b[ll] != '-' && b[ll] != '.')
                    cur_seq[h++] = lwrchr(b[ll]);
            for (; ll <= l; ll++)
                cur_seq[h++] = uprchr(b[ll]);
        }
        iprev = i;
        lprev = l;

        if (h >= (size_t)(par_maxcol - 1000))  // too few columns? Reserve double space
        {
            char *new_seq = new char[2 * par_maxcol];
            strncpy(new_seq, cur_seq, h);  //////// check: maxcol-1 ????
            delete[] (cur_seq);
            cur_seq = new_seq;
            par_maxcol *= 2;
        }
    }

    // Count match states in alignment
    int matches = 0,
        inserts = 0,
        deletions = 0;
    for (i = 0; i < result.backtrace.length(); i++) {
        switch (result.backtrace[i]) {
            case 'M':
                matches++;
                break;
            case 'I':
                inserts++;
                break;
            case 'D':
                deletions++;
                break;
        }
    }
    /* size_t numMatches = std::count(result.backtrace.begin(), result.backtrace.end(), 'M'); */
    /* size_t nonMatches = result.backtrace.length() - numMatches; */

    // Add the remaining gaps '-' to the end of the template sequence
    if (result.qEndPos > result.dbEndPos) {
        size_t remainingGaps = result.qEndPos - result.dbEndPos + inserts - deletions;
        std::cout << "Remaining: " << remainingGaps << std::endl;
        for (i = 0; i < remainingGaps; i++)
            cur_seq[h++] = '-';
        printCurSeq(cur_seq, h);
    }

    cur_seq[h++] = '\0';

    // Print out cur_seq
    std::cout
        << "               "
        << std::string(result.dbStartPos + result.qStartPos, ' ') << '*'
        << std::string(result.backtrace.length() - 2, ' ') << '*'
        << std::endl;
    std::cout << "    a3m query: " << a3m_a;
    std::cout << "\n   a3m target: " << a3m_b;
    std::cout
        << "\n               "
        << std::string(result.qStartPos + std::max(result.dbStartPos, 0), ' ') << '*'
        << std::string(result.backtrace.length() - 2, ' ') << '*'
        << std::endl;
    std::cout << "\n        Query: " << a3m_a;
    std::cout << "\n      cur_seq: ";
    for (i = 0; i < h; i++)
        std::cout << cur_seq[i];

    delete[] cur_seq;
    std::cout << std::endl;
}

int generatetree(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    // Query DB, AA and structure
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbrAA(par.db1, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader qdbr(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);

    // Alignment database from structurealign
    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    
    size_t dbSize = reader.getSize();
    unsigned int thread_idx = 0;

    std::vector<Matcher::result_t> result;
    std::vector<std::pair<size_t, Matcher::result_t > > results;    // Query ID, best hit result_t

    std::cout << "Generating tree" << std::endl;

    MSASequence sequences[dbSize];
    MSANode nodes[dbSize];

    // Save best hit per query, reverse sorted by score
    for (size_t i = 0; i < dbSize; ++i) {
        char *data = reader.getData(i, thread_idx);
        Matcher::readAlignmentResults(result, data, 0);

        // Save best hit, ignoring self-hits
        if (result.size() > 1)
            results.emplace_back(reader.getDbKey(i), result[1]);

        // Get sequence of query node
        unsigned int queryId = qdbr.sequenceReader->getId(reader.getDbKey(i));

        // FIXME
        // this is likely going out of scope, so address will be garbage when
        // stored in members array
        const char *cseq = qdbrAA.sequenceReader->getData(queryId, thread_idx);
        size_t seqLen = qdbrAA.sequenceReader->getSeqLen(queryId);
        sequences[i].id = queryId;
        sequences[i].sequence = std::string(cseq, seqLen);

        nodes[i].id = queryId;
        nodes[i].members.push_back(&sequences[i]);

        result.clear();
    }
    std::sort(results.begin(), results.end(), BestResult);

    // Find linkages
    // Does a single pass through the best hits, saving unique pairs of query and target
    std::set<int> merged;
    std::map<int, MSANode> nodies;
    std::vector<int> merges;

    //  Sequences = []
    //  Hits = []
    //  For each entry in align db (query sequence)
    //      new MSASequence (query ID, sequence from sequenceDB) --> Sequences[i]
    //      save best hit (target, dbKey) --> (MSASequence, result_t) --> Hits[i]
    //
    //  Reverse sort Hits (best --> worst)
    //
    //  Merged = {}
    //  Nodes = []
    //  For each best hit (MSASequence, result_t)
    //
    //      query ID  = MSASequence.id
    //      target ID = result_t.dbKey
    //
    //      If query ID in Merged
    //          get MSANode
    //      Else
    //          create MSANode --> Nodes
    //          MSANode --> Merged[query ID]
    //
    //      If target ID in Merged
    //          get MSANode
    //      Else
    //          create MSANode --> Nodes
    //          MSANode --> Merged[target ID]
    //
    //      If queryNode == targetNode
    //          continue
    //      
    //      # always assume result_t is hit to targetNode from queryNode
    //      # (i.e. dbKey = target ID)
    //      # Need queryID corresponding to query sequence of result_t
    //      queryNode.add_cluster(queryID, targetNode, result_t)

    MSANode *query;
    MSANode *target;

    // Create array holding MSANode pointers
    // Access the most current MSANode from any sequence id
    MSANode* nodes_p[dbSize];
    for (size_t i = 0; i < dbSize; ++i)
        nodes_p[i] = &nodes[i];
    
    for (MSANode *node : nodes_p) {
        node->print(); 
        std::cout << std::endl;
    }

    std::cout << "Merging nodes_p\n";
    for (std::pair<size_t, Matcher::result_t> pair : results) {
        size_t queryId = pair.first;
        Matcher::result_t result = pair.second;

        query = nodes_p[queryId];
        target = nodes_p[result.dbKey];

        // Check if query or target have already been merged
        if (query->id == target->id)
            continue;

        std::cout << "Merging " << result.dbKey << " (" << target->id << ")" << " into " << queryId << " (" << query->id << ")" << std::endl;
        query->addNode(queryId, target, result);
               
        // overwrite MSANode at q/t index after merge
        // so all indexes point to MSANode sequence is currently in
        for (size_t i = 0; i < target->members.size(); ++i)
            nodes_p[target->members[i]->id] = query;

        // Print
        std::cout << "\n******\n";
        std::set<int> set;
        for (size_t i = 0; i < dbSize; ++i) {
            MSANode *node = nodes_p[i]; 
            if (set.find(node->id) != set.end()) {
                //std::cout << "Found " << node->id << " in set\n";
                continue;
            }

            std::cout << "\n-----------------------------\n";
            std::cout << "Node " << node->id << std::endl;
            node->print();
            std::cout << "-----------------------------\n";

            set.insert(node->id);
            //for (int id : set)
            //    std::cout << id << " ";
            std::cout << std::endl;
        }
 
    }

    for (size_t i = 0; i < dbSize; ++i) {
        std::cout << i << " = " << nodes_p[i]->id << std::endl;
    }

    // Print the MSA
    std::set<int> set;
    for (size_t i = 0; i < dbSize; ++i) {
        MSANode *node = nodes_p[i]; 
        if (set.find(node->id) != set.end()) {
            std::cout << "Found " << node->id << " in set\n";
            continue;
        }

        std::cout << "Printing node " << node->id << std::endl;
        node->print();

        set.insert(node->id);
        for (int id : set)
            std::cout << id << " ";
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
