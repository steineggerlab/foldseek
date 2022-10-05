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

struct CompareResult {
    public:
        template <class T> bool operator() (T &a, T &b) {
            return a.second.score > b.second.score;
        }
} BestResult;

struct CompareResultLen {
    public:
        template <class T> bool operator() (T &a, T &b) {
            return a.second.alnLength > b.second.alnLength;
        }
} BestResultLen;



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
    
    // DBWriter for alignment results
    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_OMIT_FILE);
    resultWriter.open();
    
    size_t dbSize = reader.getSize();
    unsigned int thread_idx = 0;
    
    std::cout << "dbSize: " << dbSize << std::endl;

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
        if (result.size() > 1) {
            // results.emplace_back(reader.getDbKey(i), result[1]);
            for (size_t j = 1; j < result.size(); ++j) {
                results.emplace_back(reader.getDbKey(i), result[j]);
            }
        }

        // Get sequence of query node
        unsigned int queryId = qdbr.sequenceReader->getId(reader.getDbKey(i));
        std::cout << "Index: " << i << ", getDbKey: " << reader.getDbKey(i) << ", queryId: " << queryId << std::endl;
        
        const char *cseq = qdbrAA.sequenceReader->getData(queryId, thread_idx);
        size_t seqLen = qdbrAA.sequenceReader->getSeqLen(queryId);
        
        sequences[i].id = queryId;
        sequences[i].sequence = std::string(cseq, seqLen);
        sequences[i].originalLen = seqLen;

        nodes[i].id = queryId;
        nodes[i].members.push_back(&sequences[i]);

        result.clear();
    }
    std::sort(results.begin(), results.end(), BestResultLen);

    std::cout << "Found " << results.size() << " hits.\n";
    std::cout << "Found " << sizeof(sequences) / sizeof(MSASequence) << " sequences.\n";

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
        unsigned int queryId = pair.first;
        Matcher::result_t result = pair.second;

        query = nodes_p[queryId];
        target = nodes_p[result.dbKey];

         // Check if query or target have already been merged
        if (query->id == target->id)
            continue;
       
        // Always align node with less members into node with more
        // TODO should store some attribute per MSANode about no. aligned columns instead of size()
        std::cout << "Best hit query: " << queryId << ", target: " << result.dbKey << std::endl;
        int queryCount = query->countColumns(0.2);
        int targetCount = target->countColumns(0.2);
        std::cout << "Target (" << target->id << ") count: " << targetCount << std::endl;;
        std::cout << "Query  (" << query->id << ") count: " << queryCount << std::endl;
        std::cout << "Target > query? " << (targetCount > queryCount) << std::endl;
        // if (target->members.size() > query->members.size()) {
        if (targetCount > queryCount) {
            std::cout << ">>>>>>>>>>>Swapping target and query\n";
            std::swap(query, target);
            std::swap(queryId, result.dbKey);
            std::swap(result.dbStartPos, result.qStartPos);
            std::swap(result.dbEndPos, result.qEndPos);
            for (size_t i = 0; i < result.backtrace.length(); ++i) {
                if (result.backtrace[i] == 'D')
                    result.backtrace[i] = 'I';
                else if (result.backtrace[i] == 'I')
                    result.backtrace[i] = 'D';
            }
        }

        std::cout << "Merging " << result.dbKey << " (" << target->id << ")" << " into " << queryId << " (" << query->id << ")" << std::endl;
        query->addNode(queryId, target, result);
               
        // overwrite MSANode at q/t index after merge
        // so all indexes point to MSANode sequence is currently in
        for (size_t i = 0; i < target->members.size(); ++i)
            nodes_p[target->members[i]->id] = query;

        // Print
        std::set<int> set;
        for (size_t i = 0; i < dbSize; ++i) {
            MSANode *node = nodes_p[i]; 
            if (set.find(node->id) != set.end()) {
                continue;
            }

            // std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
            // std::cout << "Node " << node->id << std::endl;
            // node->print();
            // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
            // std::cout << std::endl;

            set.insert(node->id);
        }
 
    }

    // Print the MSA
    std::set<int> set;
    std::vector<int> finalNodeIds;
    for (size_t i = 0; i < dbSize; ++i) {
        MSANode *node = nodes_p[i]; 
        if (set.find(node->id) != set.end()) {
            continue;
        }
        set.insert(node->id);
        finalNodeIds.push_back(node->id);
    }
    
    resultWriter.writeStart(0);
    
    std::string entry;
    for (size_t i = 0; i < finalNodeIds.size(); ++i) {
        int nodeId = finalNodeIds[i];
        MSANode *node = nodes_p[nodeId]; 

        // Write to alignmentDB
        for (MSASequence *seq : node->members) {
            entry = seq->asString();
            resultWriter.writeAdd(entry.c_str(), entry.size(), 0);
        }
    }
    
    resultWriter.writeEnd(0, 0, false, 0);
    resultWriter.close(true);
    FileUtil::remove(par.db3Index.c_str());

    return EXIT_SUCCESS;
}
