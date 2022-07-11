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
#include <tuple>
#include <set>

#ifdef OPENMP
#include <omp.h>
#endif

#define	EXIT_FAILURE	1
#define	EXIT_SUCCESS	0

struct Node {
    size_t id;
    std::string sequence;
    Matcher::result_t result;

    Node(size_t n_id, std::string n_sequence, Matcher::result_t n_result)
        : id(n_id), sequence(n_sequence), result(n_result) {};
    Node() = default;

    void print() {
        std::cout << ">Node " << id << "\n" << sequence;
    }
};

/**
 * A structure cluster within the MSA.
 * Stores member Nodes and a representative sequence for adding new Nodes.
 */
struct Cluster {
    size_t id;
    std::vector<Node *> members;
    std::string representative;
};

struct Linkage {
    size_t query;
    size_t target;
    int score;
    size_t index;
};

struct CompareResult{
    public:
        template <class T> bool operator() (T &a, T &b) {
            return a.second.score > b.second.score;
        }
} BestResult;


template <class T>
void printVector(std::vector<T> vector) {
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
    /* std::cout << "  " << a; */
    /* std::cout << "  " << b; */
}

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


void testing(std::string a, std::string b, Matcher::result_t result) {

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

    int par_maxcol = a.length() + b.length();
    char *cur_seq = new char[par_maxcol];
    char c;
    size_t h;
    size_t l, ll;
    int pos;

    // Add left-end gaps
    for (h = 0; h < (size_t)result.qStartPos; h++)
        cur_seq[h] = '-';

    // Mark insertions until first match position
    j = 0;
    pos = 0;
    if (result.dbStartPos > 0) {
        while (j < (size_t)result.dbStartPos && b[pos] != '\0') {
            c = b[pos++];
            if (match(c)) j++;
            cur_seq[h++] = (c == '-') ? '.' : lwrchr(c);
        }
    }

    std::cout << "Marked insertions:\n";
    i = 0;
    while (cur_seq[i] != '\0') {
        std::cout << cur_seq[i];
        i++;
    }
    std::cout << "\n\n\n";

    // Check if start position is an insertion; if so, copy
    // b[pos] is first position after marking insertions above
    c = b[pos];
    while (insertion(c) && (c != '\0')) {
        cur_seq[h++] = lwrchr(c);
        c = b[pos++];
    }

    std::cout << "Checked start position:\n";
    i = 0;
    while (cur_seq[i] != '\0') {
        std::cout << cur_seq[i];
        i++;
    }
    std::cout << "\n\n\n";


    // Advance to match state result.dbStartPos of b
    // match state at position l?
    // yes: increment j. Reached hit,j1? yes: break
    // TODO is this any different to just l = result.dbStartPos ?
    //  Just use pos?

    for (j = 0, l = 0; (c = b[l]) != '\0'; l++) {
        if (match(c))
            if ((++j) == (size_t)result.dbStartPos)
                break;
    }
    /* l = result.dbStartPos - 1; */


    // Write first match state to cur_seq
    int iprev = result.qStartPos;  // Previous query match state
    int lprev = l;                 // Previous target match state
    cur_seq[h++] = b[l];        // First column is Match-Match state

    // Iterate next match states
    for (j = result.dbStartPos + 1; j <= (result.dbStartPos + result.backtrace.length()); j++) {
        // Get index in query from imap
        i = imap[j];

        // Advance to position of next T match state j
        while ((c = b[++l]) > '\0' && insertion(c));

        int di = i - iprev;
        int dl = l - lprev;

        /* std::cout */
        /*     << "i: " << i << ", iprev: " << iprev << ", di: " << di */
        /*     << ", l: " << l << ", lprev: " << lprev << ", dl: " << dl */
        /*     << std::endl; */
        
        if (di == 1) {
            /* std::cout << "di == 1" << std::endl; */
            for (ll = lprev + 1; ll < l; ll++)
                if (b[ll] != '-' && b[ll] != '.')
                    cur_seq[h++] = lwrchr(b[ll]);
            cur_seq[h++] = b[ll++];

        } else if (di == 0) {
            /* std::cout << "di == 0" << std::endl; */
            for (ll = lprev + 1; ll <= l; ll++)
                if (b[ll] != '-' && b[ll] != '.')
                    cur_seq[h++] = lwrchr(b[ll]);

        } else if (di >= dl) {
            /* std::cout << "di >= dl" << std::endl; */
            for (ll = lprev + 1; ll <= (size_t)(lprev + (int)(dl / 2)); ll++)
                cur_seq[h++] = uprchr(b[ll]);
            for (int gap = 1; gap <= di - dl; gap++)
                cur_seq[h++] = '-';
            for (; ll <= l; ll++)
                cur_seq[h++] = uprchr(b[ll]);

        } else if (di < dl) {
            /* std::cout << "di < dl" << std::endl; */
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

        if (h >= par_maxcol - 1000)  // too few columns? Reserve double space
        {
            char *new_seq = new char[2 * par_maxcol];
            strncpy(new_seq, cur_seq, h);  //////// check: maxcol-1 ????
            delete[] (cur_seq);
            cur_seq = new_seq;
            par_maxcol *= 2;
        }

        i = 0;
        while (cur_seq[i] != '\0') {
            std::cout << cur_seq[i];
            i++;
        }
        std::cout << '\n';
    }

    // Count match states in alignment
    size_t L = std::count(result.backtrace.begin(), result.backtrace.end(), 'M');

    // Add the remaining gaps '-' to the end of the template sequence
    for (i = result.qEndPos + 1; i <= L; i++) {
        cur_seq[h++] = '-';

    }

    // add remaining seq. info as insertion state
    while (b[ll] != '\0') {
        if (b[ll] == '-') {
            cur_seq[h++] = '.';
        } else {
            cur_seq[h++] = lwrchr(a[ll]);
        }
        ll++;
    }

    cur_seq[h++] = '\0';


    // Print out cur_seq
    {
        std::cout
            << "               "
            << std::string(result.dbStartPos + result.qStartPos, ' ') << '*'
            << std::string(result.backtrace.length() - 2, ' ') << '*'
            << std::endl;
        std::cout << "    a3m query: " << a3m_a;
        std::cout << "   a3m target: " << a3m_b;
        std::cout
            << "               "
            << std::string(result.qStartPos + std::max(result.dbStartPos, 0), ' ') << '*'
            << std::string(result.backtrace.length() - 2, ' ') << '*'
            << std::endl;
        std::cout << "        Query: " << a3m_a;
        std::cout << "      cur_seq: "; //<< std::string(result.dbStartPos, ' ');
        i = 0;
        while (cur_seq[i] != '\0') {
            std::cout << cur_seq[i];
            i++;
        }
    }

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

    std::vector<int> id2cluster(dbSize, -1);  // Map database IDs to cluster IDs
    std::vector<Node> nodes;

    // Save best hit per query, reverse sorted by score
    for (size_t i = 0; i < dbSize; ++i) {
        char *data = reader.getData(i, thread_idx);
        Matcher::readAlignmentResults(result, data, 0);

        // Save best hit, ignoring self-hits
        if (result.size() > 1)
            results.emplace_back(reader.getDbKey(i), result[1]);

        std::cout << "Align length: " << result[1].alnLength << ", Backtrace: " << result[1].backtrace.length() << std::endl;
        
        // Get sequence of query node
        unsigned int queryId = qdbr.sequenceReader->getId(reader.getDbKey(i));
        nodes.emplace_back(reader.getDbKey(i), qdbrAA.sequenceReader->getData(queryId, thread_idx), result[1]);

        result.clear();
    }
    std::sort(results.begin(), results.end(), BestResult);

    // TODO: make a simpler test case for a3m merging
    // e.g.
    //  Q   X X X X X X X X X X X X X X X
    //  T1  Y Y Y Y Y Y Y Y Y Y Y
    //  T2  Z Z Z Z Z Z Z Z
    //
    //  Q-T1 CIGAR
    //        M M M D D M M M M M M M M
    //  Q   X X X X X X X X X X X X X X X
    //  T1    Y Y Y - - Y Y Y Y Y Y Y Y
    //
    //  Q-T2
    //          M M M M M I I M
    //  T1  Y Y Y Y Y Y Y . . Y Y Y Y
    //  T2      Z Z Z Z Z z z Z
    //
    //  Expected a3m:
    //  Q   X X X X X X X X X X X X X X X
    //  T1  - Y Y Y - - Y Y Y Y Y Y Y Y
    //  T2  - - - Z - - Z Z Z Z z z Z
    //
    //  Full alignment:
    //  Q   X X X X X X X X X X - - X X X X X
    //  T1  - Y Y Y - - Y Y Y Y - - Y Y Y Y -
    //  T2  - - - Z - - Z Z Z Z z z Z - - - -
    //
    //  TODO: use CIGAR strings from other sequences to inform exact positioning ?

    /* unsigned int dbKey; */
    /* int score; */
    /* float qcov; */
    /* float dbcov; */
    /* float seqId; */
    /* double eval; */
    /* unsigned int alnLength; */
    /* int qStartPos; */
    /* int qEndPos; */
    /* unsigned int qLen; */
    /* int dbStartPos; */
    /* int dbEndPos; */
    /* unsigned int dbLen; */
    /* int queryOrfStartPos; */
    /* int queryOrfEndPos; */
    /* int dbOrfStartPos; */
    /* int dbOrfEndPos; */
    /* std::string backtrace; */

    /* nodes.emplace_back(0, "XXXXXXXXXXXXXXX\n", *(new Matcher::result_t())); */
    /* nodes.emplace_back(1, "YYYYYYYYYYY\n", *(new Matcher::result_t())); */
    /* nodes.emplace_back(2, "ZZZZZZZZ\n", *(new Matcher::result_t())); */

    /* nodes[0].result.dbKey = 1; */
    /* nodes[0].result.backtrace = "MMMIIMMMMMMMM"; */
    /* nodes[0].result.qStartPos = 1; */
    /* nodes[0].result.qEndPos = 14; */
    /* nodes[0].result.dbStartPos = 0; */
    /* nodes[0].result.dbEndPos = 11; */

    /* nodes[1].result.dbKey = 2; */
    /* nodes[1].result.backtrace = "MMMMMDDM"; */
    /* nodes[1].result.qStartPos = 2; */
    /* nodes[1].result.qEndPos = 7; */
    /* nodes[1].result.dbStartPos = 0; */
    /* nodes[1].result.dbEndPos = 7; */

    /* merge_a3m(nodes[0].sequence, nodes[1].sequence, nodes[0].result); */

    std::cout << "Best hits:" << std::endl;
    for (size_t i = 0; i < results.size(); i++) {
        std::cout << "Merging " << results[i].first << " vs " << results[i].second.dbKey << std::endl;
        merge_a3m(nodes[results[i].first].sequence, nodes[results[i].second.dbKey].sequence, results[i].second);
        std::cout << std::endl;
        /* break; */
    }
    
    // Find linkages
    // Does a single pass through the best hits, saving unique pairs of query and target
    std::vector<Linkage> merges;

    /* size_t idx = results.size();  // TODO need to get maximum ID in database+1 */
    /* for (size_t i = 0; i < results.size(); i++) { */
    /*     if (merged.find(i) != merged.end() || merged.find(results[i].dbKey) != merged.end()) { */
    /*         continue; */
    /*     } else { */
    /*         Linkage linkage = {i, results[i].dbKey, results[i].score, idx}; */

    /*         merges.push_back(linkage); */
    /*         merged.insert(i); */
    /*         merged.insert(results[i].dbKey); */
    /*         idx++; */
    /*     } */
    /* } */

    /* std::cout << "\nMerged nodes:" << std::endl; */
    /* for (auto it = merged.begin(); it != merged.end(); it++) { */
    /*     std::cout << "  " << *it; */
    /* } */

    /* std::cout << "\nMerges:" << std::endl; */
    /* for (size_t i = 0; i < merges.size(); i++) { */
    /*     std::cout << "  " << merges[i].query << " " << merges[i].target << " " << merges[i].score << " " << merges[i].index << std::endl; */
    /* } */

    return EXIT_SUCCESS;
}
