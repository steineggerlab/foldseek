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
        std::cout << ">Node " << id << "\n" << sequence << std::endl;
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

void print_sequences(std::string a, std::string b) {
    std::cout << "  " << a;
    std::cout << "  " << b;
}

// Tests c is a match
inline bool match(char c) {
    return (c >= 'A' && c <= 'Z') || c == '-';
}

// Transforms chr into an uppercase character
inline char uprchr(char chr) {
  return (chr >= 'a' && chr <= 'z') ? chr + 'A' - 'a' : chr;
}

// Transforms chr into lowercase character
inline char lwrchr(char chr) {
    return (chr >= 'A' && chr <= 'Z') ? chr - 'A' + 'a' : chr;
}

/**
 * Merge alignment B into alignment A.
 */
void merge_a3m(std::string a, std::string b, Matcher::result_t result) {
    std::cout << "Merging:" << std::endl;
    std::cout << "  qStart/qEnd: " << result.qStartPos << " " << result.qEndPos << std::endl;
    std::cout << "  tStart/tEnd: " << result.dbStartPos << " " << result.dbEndPos << std::endl;
    print_sequences(a, b);

    std::cout << std::endl;

    int j;


    // Count match states (capitals) in alignment A
    // Used when adding remaining gaps to end of B
    int l, L;
    for (L = 0, l = 1; a[l] > '\0'; l++) {
        if (match(a[l])) L++;
    }

    // For each sequence in B, align to A
    int par_maxcol = a.length();
    char *cur_seq = new char[par_maxcol];

    {
        // NOTE: Treat qStartPos-qEndPos/dbStartPos-dbEndPos as match states
        
        // Lower case characters in A, if necessary
        int i;
        for (i = 0; i < result.qStartPos - 1; i++)
            a[i] = lwrchr(a[i]);
        for (i = 0; i < result.dbStartPos - 1; i++)
            b[i] = lwrchr(b[i]);
        for (i = result.qEndPos + 1; i < a.length(); i++)
            a[i] = lwrchr(a[i]);
        for (i = result.dbEndPos + 1; i < b.length(); i++)
            b[i] = lwrchr(b[i]);

        // Insert 5' gaps
        if (result.qStartPos == result.dbStartPos) {
        } else if (result.qStartPos > result.dbStartPos) {
            b.insert(0, result.qStartPos - result.dbStartPos, '-');
        } else {
            a.insert(0, result.dbStartPos - result.qStartPos, '.');
        }


        cur_seq[0] = ' ';  // 0th position not used
        
        // Add qStartpos - 1 left-end gaps to aligned sequence
        int h;
        for (h = 1; h < result.dbStartPos; h++)
            cur_seq[h] = '-';

        // 
        char c;
        j = 1;
        int pos = 1; 
        if (result.dbStartPos > 1) {
            while (j < result.dbStartPos && a[pos] != '\0') {
                c = a[pos];
                if (match(c)) j++;                    
                cur_seq[h++] = (c == '-') ? '.' : lwrchr(c);
                pos++;
            }
        }

        // Check if start position is an insertion
        // If yes, copy
        c = b[pos];
        while (((c >= 'a' && c <= 'z') || c == '.') && c != '\0') {
            cur_seq[h++] = lwrchr(c);
            pos++;
            c = b[pos];
        }

        // Advance to match state dbStartPos of b
        for (j = 0, l = 1; (c = b[l]) > '\0'; l++) {
            if (match(c)) {
                // Match at position l? increment j
                // Reached dbStartPos? break
                if ((++j) == result.dbStartPos) break;
            }
        }

        if (j < result.dbStartPos) {
            std::cerr << "Error: did not find " << result.dbStartPos << " match states";
            exit(1);
        }

        /* // Write first match state to cur_seq */
        /* int iprev = hit.i1;  // index of previous query match state */
        /* int lprev = l;      // previous T match state in Tali.seq[k][l] */
        /* cur_seq[h++] = Tali.seq[k][l];  // first column of alignment is Match-Match state */

        /* // For each further match state j in alignment */
        /* step = hit.nsteps; */
        /* for (j = hit.j1 + 1; j <= hit.j2; ++j) { */
        /*     // Advance to position of next T match state j */
        /*     i = imatch[j]; */

        /*     // Advance to position of next T match state j */
        /*     while ((c = Tali.seq[k][++l]) > '\0' */
        /*            && ((c >= 'a' && c <= 'z') || c == '.')); */

        /*     int di = i - iprev;  // number of Match states in Q between T match state j-1 and j */
        /*     int dl = l - lprev;  // 1 + number of inserted residues in T sequence between T match state j-1 and j */
        /*     if (di == 1) { */
        /*         // One Q match state for one T match state (treated as special case for speed reasons) */
        /*         // i:       i-1   i         di=1 */
        /*         // Q:  XXXXXX.....XXXXXX */
        /*         // T:  YYYYYYyyyyyYYYYYY */
        /*         // j:       j-1   j */
        /*         // l:       lprev l         dl=6 */

        /*         // Inserts in lower case */
        /*         for (ll = lprev + 1; ll < l; ll++) */
        /*             if (Tali.seq[k][ll] != '-' && Tali.seq[k][ll] != '.') */
        /*                 cur_seq[h++] = lwrchr(Tali.seq[k][ll]); */

        /*         // Template Match state -> upper case */
        /*         cur_seq[h++] = Tali.seq[k][ll]; */
        /*         ll++; */
        /*     } else if (di == 0) { */
        /*         // Gap in query: no Q match state for on T match state (special case for speed reasons) */
        /*         // i:       i-1   i-1       di=0 */
        /*         // Q:  XXXXXX.....~~~XXX */
        /*         // T:  YYYYYYyyyyyYYYYYY */
        /*         // j:       j-1   j */
        /*         // l:       lprev l         dl=6 */

        /*         // All T residues (including T match state) in lower case */
        /*         for (ll = lprev + 1; ll <= l; ll++) */
        /*             if (Tali.seq[k][ll] != '-' && Tali.seq[k][ll] != '.') */
        /*                 cur_seq[h++] = lwrchr(Tali.seq[k][ll]); */
        /*     } else if (di >= dl) { */
        /*         // More Match states in Q than Inserts in the T sequence */
        /*         // => half T inserts y left, half right-aligned in uc, gaps to fill up */
        /*         // Number of T insert residues to be left-aligned: (int)(dl/2) */
        /*         // i:        iprev  i       di=7 */
        /*         // Q:  XXXXXXXXXXXXXXXXXX */
        /*         // T:  YYYYYYYyyy-yyYYYYY */
        /*         // j:        j-1    j */
        /*         // l:        lprev  l       dl=6 */

        /*         // Add left-bounded template residues */
        /*         for (ll = lprev + 1; ll <= lprev + (int) (dl / 2); ll++) */
        /*             cur_seq[h++] = uprchr(Tali.seq[k][ll]); */

        /*         // Add central gaps */
        /*         for (int gap = 1; gap <= di - dl; gap++) */
        /*             cur_seq[h++] = '-'; */

        /*         // Add right-bounded residues */
        /*         for (; ll <= l; ll++) */
        /*             cur_seq[h++] = uprchr(Tali.seq[k][ll]); */
        /*     } else if (di < dl) { */
        /*         // Fewer Match states in Q than inserts in T sequence */
        /*         // => half of available space di for left- half for right-aligned T inserts, rest in lc */
        /*         // number of T inserts to be left-aligned in uc: (int)(di/2), */
        /*         // i:        iprev i       di=5 */
        /*         // Q:  XXXXXXXXX.XXXXXXX */
        /*         // T:  YYYYYYYyyyyyYYYYY */
        /*         // j:        j-1   j */
        /*         // l:        lprev l       dl=6 */

        /*         // Add left-bounded template residues */
        /*         for (ll = lprev + 1; ll <= lprev + (int) (di / 2); ll++) */
        /*             cur_seq[h++] = uprchr(Tali.seq[k][ll]); */

        /*         // Add central inserts */
        /*         for (int ins = 1; ins <= dl - di; ins++, ll++) */
        /*             if (Tali.seq[k][ll] != '-' && Tali.seq[k][ll] != '.') */
        /*                 cur_seq[h++] = lwrchr(Tali.seq[k][ll]); */

        /*         // Add right-bounded residues */
        /*         for (; ll <= l; ll++) */
        /*             cur_seq[h++] = uprchr(Tali.seq[k][ll]); */
        /*     } */
        /*     HH_LOG(DEBUG3) << "i=" << i << " j=" << j << " l=" << l << " cur_seq=" */
        /*                    << cur_seq << "\n"; */

        /*     iprev = i; */
        /*     lprev = l; */
        /*     if (h >= maxcol - 1000)  // too few columns? Reserve double space */
        /*     { */
        /*         char *new_seq = new char[2 * maxcol]; */
        /*         strncpy(new_seq, cur_seq, h);  //////// check: maxcol-1 ???? */
        /*         delete[] (cur_seq); */
        /*         cur_seq = new_seq; */
        /*         maxcol *= 2; */
        /*     } */
        /* } */


        print_sequences(a, b);

    }


    /* print_sequences(a, b); */

    /* int i = 1; */
    /* int pos = 1; */
    /* char c; */
    /* while (i < result.qStartPos && b[pos] != '\0') { */
    /*     std::cout << i << " " << b[pos] << " "; */
    /*     c = b[pos]; */

    /*     if ((c >= 'A' && c <= 'Z') || c == '-') { */
    /*         i++; */
    /*     } */

    /*     b[pos++] = */ 

    /*     pos++; */
    /* } */

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
        
        // Get sequence of query node
        unsigned int queryId = qdbr.sequenceReader->getId(reader.getDbKey(i));
        nodes.emplace_back(reader.getDbKey(i), qdbrAA.sequenceReader->getData(queryId, thread_idx), result[1]);

        result.clear();
    }
    std::sort(results.begin(), results.end(), BestResult);

    for (size_t i = 0; i < nodes.size(); i++) {
        nodes[i].print();
    }

    std::cout << "Best hits:" << std::endl;
    for (size_t i = 0; i < results.size(); i++) {
        std::cout << results[i].first << " vs " << results[i].second.dbKey << " " << results[i].second.score << std::endl;

        merge_a3m(nodes[results[i].first].sequence, nodes[results[i].second.dbKey].sequence, results[i].second);
    }
    
    size_t maxId = dbSize + 1;
    std::set<int> merged;
    std::vector<Cluster> clusters;
    for (size_t i = 0; i < results.size(); i++) {
        bool qMerged = (id2cluster[results[i].first] != results[i].first);
        bool tMerged = (id2cluster[results[i].second.dbKey] != results[i].second.dbKey);
        if (!qMerged && !tMerged) {
            std::vector<size_t> members{results[i].first, results[i].second.dbKey};
            /* clusters.emplace_back(maxId, members, ); */
            maxId++;
        } else {
        }
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
