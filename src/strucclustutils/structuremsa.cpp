#include "Alignment.h"
#include "BacktraceTranslator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "IndexReader.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "MathUtil.h"
#include "MsaFilter.h"
#include "MultipleAlignment.h"
#include "Newick.h"
#include "PSSMCalculator.h"
#include "Parameters.h"
#include "Sequence.h"
#include "StructureSmithWaterman.h"
// #include "affineneedlemanwunsch.h"
#include "StructureUtil.h"
#include "Util.h"
#include "structureto3diseqdist.h"
#include <cassert>
#include <tuple>
#include <set>
#include <fstream>
#include <iostream>
#include <regex>
#include <stack>

#include "kseq.h"
#include "KSeqBufferReader.h"
#include "LDDT.h"

#ifdef OPENMP
#include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define	EXIT_FAILURE	1
#define	EXIT_SUCCESS	0

struct AlnSimple {
    unsigned int queryId;
    unsigned int targetId;
    int score;
};

struct TNode {
    int id;
    int dbId;
    std::string name;
    float length;
    std::vector<TNode *> children;
    TNode *parent;
};

void printTree(TNode *node, int depth) {
    std::string gap(depth * 2, ' ');
    std::cout << gap << node->id;
    if (node->name != "") {
        std::cout << "=" << node->name;
        std::cout << " (parent " << node->parent->id << ", length " << node->length;
        if (node->dbId != -1) {
            std::cout << ", dbId " << node->dbId;
        }
        std::cout << ")";
    }
    std::cout << '\n';
    for (size_t i = 0; i < node->children.size(); i++) {
        TNode *child = node->children[i];
        printTree(child, depth + 1);
    }
}

/**
 * @brief Post-order traversal of a parsed Tree.
 * Generates the merging order for structuremsa
 * 
 * @param node Pointer to root TNode of the tree
 */
void postOrder(TNode *node, std::vector<int> *linkage) {
    for (TNode *child : node->children) {
        postOrder(child, linkage);
    }
    if (node->children.size() > 0) {
        for (TNode *child : node->children) {
            linkage->push_back(child->dbId);

            // Propagate child dbId from leaf to root, so we
            // always have a reference during alignment stage
            node->dbId = child->dbId;
        }
    }
}

std::vector<AlnSimple> parseNewick(std::string newick, std::map<std::string, int> &headers) {
    // Should know this from number of structures in database
    int total = std::count(newick.begin(), newick.end(), '(');

    // Tokenize tree on ; | ( ) , :
    // Use {-1, 0} on sregex_token_iterator to get matches AND inbetween (i.e. names)
    std::regex pattern("\\s*(;|\\(|\\)|,|:)\\s*");
    auto words_begin = std::sregex_token_iterator(newick.begin(), newick.end(), pattern, {-1, 0});
    auto words_end = std::sregex_token_iterator();
    
    // Initialise TNode array (2n+1 but a +1 for the root)
    TNode nodes[total * 2 + 2];
    std::vector<TNode *> parents;
    TNode *tree;
    TNode *subtree;
    
    // Initialise the root node (the +1)
    TNode root;
    root.id = 0;
    nodes[0] = root;
    tree = &nodes[0];
    
    int count = 1;
    std::string prevToken;
    
    for (std::sregex_token_iterator i = words_begin; i != words_end; ++i) {
        std::string match_str = *i;

        if (match_str == "(") {
            // add new node, set it as new subtree
            TNode newNode;
            newNode.id = count;
            newNode.dbId = -1;
            newNode.length = 1;
            nodes[count] = newNode;
            subtree = &nodes[count];
            count++;
            
            // add it as child to current tree
            tree->children.push_back(subtree);
            subtree->parent = tree;
            
            // add the tree as parent, set subtree as tree
            parents.push_back(tree);
            tree = subtree;
        } else if (match_str == ",") {
            TNode newNode;
            newNode.id = count;
            newNode.dbId = -1;
            newNode.length = 1;
            nodes[count] = newNode;
            subtree = &nodes[count];
            count++;
            parents.back()->children.push_back(subtree);
            subtree->parent = parents.back();
            tree = subtree;
        } else if (match_str == ")") {
            tree = parents[parents.size() - 1];
            parents.pop_back();
        } else if (match_str == ":") {
            // Don't do anything here, just catch it in case there are
            // branch lengths to set in else
        } else {
            if (match_str != "" && (prevToken == ")" || prevToken == "(" || prevToken == ",")) {
                tree->name = match_str;
                tree->dbId = headers[match_str];
            } else if (prevToken == ":") {
                tree->length = std::stof(match_str);
            }
        }
        prevToken = match_str;
    }
   
    // printTree(tree, 0);
    
    // Get (flat) linkage matrix, 2(2n+1)
    // node 1, node 2
    // NOTE: postOrder will trip up when no. children != 2
    //       will get duplicate rows which cause errors
    std::vector<int> linkage;
    std::vector<AlnSimple> hits;
    postOrder(tree, &linkage);
    for (size_t i = 0; i < linkage.size(); i += 2) {
        AlnSimple hit;
        hit.queryId = linkage[i + 0];
        hit.targetId = linkage[i + 1];
        hits.push_back(hit);
    }

    return hits;
}

Matcher::result_t pairwiseAlignment(StructureSmithWaterman & aligner, unsigned int querySeqLen,  Sequence *target_aa, Sequence *target_3di, int gapOpen,
                  int gapExtend) {
    std::string backtrace;
    
    bool targetIsProfile = (Parameters::isEqualDbtype(target_aa->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));

    unsigned char * target_aa_seq = target_aa->numSequence;
    unsigned char * target_3di_seq = target_3di->numSequence;
    
    if (targetIsProfile) {
        target_aa_seq = target_aa->numConsensusSequence;
        target_3di_seq = target_3di->numConsensusSequence;
    }
    StructureSmithWaterman::s_align align = aligner.alignScoreEndPos(
        target_aa_seq,
        target_3di_seq,
        target_aa->L,
        gapOpen,
        gapExtend,
        querySeqLen / 2
    );
    align = aligner.alignStartPosBacktrace(
        target_aa_seq,
        target_3di_seq,
        target_aa->L,
        gapOpen,
        gapExtend,
        3,
        backtrace,
        align,
        0,
        0.0,
        querySeqLen / 2
    );
    unsigned int alnLength = Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1, align.dbStartPos1, align.dbEndPos1);
    alnLength = backtrace.size();
    float seqId = Util::computeSeqId(Parameters::SEQ_ID_ALN_LEN, align.identicalAACnt, querySeqLen, target_aa->L, alnLength);
    return Matcher::result_t(
        target_aa->getDbKey(),
        align.score1,
        align.qCov,
        align.tCov,
        seqId,
        align.evalue,
        alnLength,
        align.qStartPos1,
        align.qEndPos1,
        querySeqLen,
        align.dbStartPos1,
        align.dbEndPos1,
        target_aa->L,
        backtrace
    );
}

void sortHitsByScore(std::vector<AlnSimple> &hits) {
    std::stable_sort(hits.begin(), hits.end(), [](const AlnSimple & a, const AlnSimple & b) {
        return a.score > b.score;
    });
}

std::vector<AlnSimple> removeMergedHits(std::vector<AlnSimple> & hits, unsigned int mergedId, unsigned int targetId) {
    std::vector<AlnSimple> newHits;
    for (size_t i = 0; i < hits.size(); i++) {
        if (hits[i].queryId != mergedId && hits[i].targetId != mergedId
            && hits[i].queryId != targetId && hits[i].targetId != targetId) {
            newHits.push_back(hits[i]);
        }
    }
    return newHits;
}


inline int get1dIndex(int i, int j, int N) {
    return j + i * (2 * N - i - 1) / 2 - i - 1;
}

std::vector<AlnSimple> updateAllScores(
    int8_t * tinySubMatAA,
    int8_t * tinySubMat3Di,
    SubstitutionMatrix * subMat_aa,
    std::vector<Sequence*> & allSeqs_aa,
    std::vector<Sequence*> & allSeqs_3di,
    bool * alreadyMerged,
    int maxSeqLen,
    int alphabetSize,
    int compBiasCorrection,
    int compBiasCorrectionScale
) {
    unsigned int N = allSeqs_aa.size();
    std::vector<AlnSimple> newHits((N * N - N) / 2);
#pragma omp parallel
{
    StructureSmithWaterman structureSmithWaterman(maxSeqLen, alphabetSize, compBiasCorrection, compBiasCorrectionScale);
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < N; i++) {
        if (alreadyMerged[i])
            continue;
        structureSmithWaterman.ssw_init(
            allSeqs_aa[i],
            allSeqs_3di[i],
            tinySubMatAA,
            tinySubMat3Di,
            subMat_aa
        );
        for (unsigned int j = i + 1; j < N; j++) {
            if (alreadyMerged[j] || i == j)
                continue;
            bool targetIsProfile = (Parameters::isEqualDbtype(allSeqs_aa[j]->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));
            unsigned char *target_aa_seq = allSeqs_aa[j]->numSequence;
            unsigned char *target_3di_seq = allSeqs_3di[j]->numSequence;
            if (targetIsProfile) {
                target_aa_seq = allSeqs_aa[j]->numConsensusSequence;
                target_3di_seq = allSeqs_3di[j]->numConsensusSequence;
            }
            AlnSimple aln;
            aln.queryId = i;
            aln.targetId = j;
            aln.score = structureSmithWaterman.ungapped_alignment(target_aa_seq, target_3di_seq, allSeqs_aa[j]->L);
            int ij = get1dIndex(i, j, N);
            newHits[ij] = aln;
        }
    }
}
    return newHits;
}

int findRoot(int vertex, std::vector<int>& parent) {
    while (parent[vertex] != vertex) {
        parent[vertex] = parent[parent[vertex]];
        vertex = parent[vertex];
    }
    return vertex;
}

TNode* findRoot(TNode *node) {
    TNode* root = node;
    while (root->id != root->parent->id)
        root = root->parent;
    return root;
}

std::string getNewick(TNode* node) {
    if (node->children.empty()) {
        return node->name;
    } else {
        std::string s = "(";
        for (size_t i = 0; i < node->children.size(); i++) {
            if (i > 0)
                s += ",";
            s += getNewick(node->children[i]);
        }
        s += ")";
        return s;
    }
}

/**
 * @brief Generates a Newick format string from a linkage matrix.
 * 
 * @param hits linkage matrix
 * @param headers structure headers
 * @param n number of structures
 * @return std::string 
 */
std::string orderToTree(std::vector<AlnSimple> hits, std::string *headers, int n) {
    std::string nwk = "";
    int totalNodes = n * 2 + 2;
    TNode nodes[totalNodes];
    TNode *root = &nodes[0];
    for (int i = 0; i < totalNodes; i++) {
        nodes[i].id = i;
        nodes[i].parent = &nodes[i];
        if (i < n) {
            nodes[i].name = headers[i];
        }
    }
    int newId = n + 1;
    for (AlnSimple aln : hits) {
        TNode *u = findRoot(&nodes[aln.queryId]);
        TNode *v = findRoot(&nodes[aln.targetId]);
        TNode *newNode = &nodes[newId];
        newNode->id = newId;
        newNode->length = aln.score;
        newNode->children.push_back(u);            
        newNode->children.push_back(v);            
        u->parent = newNode;
        v->parent = newNode;
        newId++;
        root = newNode;
    }
    // printTree(root, 0);
    return getNewick(root);
}

/**
 * @brief Get minimum spanning tree as linkage matrix (Kruskal algorithm).
 * 
 * @param hits all hits from UpdateAllScores
 * @param n number of structures
 * @return std::vector<AlnSimple> 
 */
std::vector<AlnSimple> mst(std::vector<AlnSimple> hits, int n) {
    std::vector<AlnSimple> result;
    std::vector<int> parent(n);  // parent node IDs
    std::vector<int> counts(n);  // how many times we have seen this id
    for (int i = 0; i < n; i++)
        parent[i] = i;
    for (AlnSimple aln : hits) {
        int u = findRoot(aln.queryId, parent);
        int v = findRoot(aln.targetId, parent);
        if (u != v) {
            result.push_back(aln);
            parent[u] = v;
        }
    }
    return result;
}

/**
 * @brief Reorder linkage matrix to maximize unique merges per iteration for multithreading.
 * 
 * @param linkage linkage matrix generated by `mst`
 * @param merges number of unique merges per iteration
 * @param n number of structures
 * @return std::vector<AlnSimple> 
 */
std::vector<AlnSimple> reorderLinkage(std::vector<AlnSimple> linkage, std::vector<int> &merges, int n) {
    std::vector<int> parent(n); 
    std::vector<int> counts(n);
    for (int i = 0; i < n; i++) {
        parent[i] = i;
        counts[i] = 0;
    }
    std::vector<AlnSimple> result(linkage.size());
    std::vector<bool> merged(linkage.size());
    int index = 0;
    int mergeCount = 0; // no. total merges
    int mergeTally = 0; // count per round
    while (mergeCount < (int)linkage.size()) {
        for (int i = 0; i < n; i++)
            counts[i] = 0;
        for (size_t i = 0; i < linkage.size(); i++) {
            if (merged[i])
                continue;
            AlnSimple aln = linkage[i];
            int u = findRoot(aln.queryId, parent);
            int v = findRoot(aln.targetId, parent);
            if (counts[u] > 0 || counts[v] > 0)
                continue;
            result[index++] = aln;
            parent[u] = v;
            merged[i] = true;
            counts[u]++;
            counts[v]++;
            mergeTally++;
        }
        merges.push_back(mergeTally);
        mergeCount += mergeTally;
        mergeTally = 0;
    }
    return result;
}


std::string fastamsa2profile(std::string & msa, PSSMCalculator &pssmCalculator, MsaFilter &filter, SubstitutionMatrix &subMat, size_t maxSeqLength, size_t maxSetSize,
                             float matchRatio, bool filterMsa, bool compBiasCorrection, std::string & qid, float filterMaxSeqId, float Ndiff, float covMSAThr,
                             float qsc, int filterMinEnable, bool wg, bool *externalMaskedColumns, float scoreBias) {
    enum {
        MSA_CA3M = 0,
        MSA_A3M  = 1,
        MSA_STOCKHOLM = 2
    };
    // set up parser
    kseq_buffer_t d;
    d.buffer = (char*)msa.c_str();
    d.length = msa.size();

    // filter parameter
    std::vector<std::string> qid_str_vec = Util::split(qid, ",");
    std::vector<int> qid_vec;
    for (size_t qid_idx = 0; qid_idx < qid_str_vec.size(); qid_idx++) {
        float qid_float = strtod(qid_str_vec[qid_idx].c_str(), NULL);
        qid_vec.push_back(static_cast<int>(qid_float*100));
    }
    std::sort(qid_vec.begin(), qid_vec.end());

    // default parameter
    bool fastaError = false;
    bool maskByFirst = false;
    kseq_t *seq = kseq_init(&d);
    // bool inHeader = false;
    unsigned int setSize = 0;
    // unsigned int seqLength = 0;
    size_t msaPos = 0;
    unsigned int centerLengthWithGaps = 0;
    unsigned int maskedCount = 0;
    unsigned int msaType = 2; // stockholm

    // init memory
    bool *maskedColumns = new bool[maxSeqLength + 1];
    Sequence sequence(maxSeqLength + 1, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, false);
    char **msaSequences = (char**) mem_align(ALIGN_INT, sizeof(char*) * maxSetSize);
    char *msaContent = (char*) mem_align(ALIGN_INT, sizeof(char) * (maxSeqLength + 1) * maxSetSize);
    float *seqWeight = new float[maxSetSize];
    float *pNullBuffer = new float[maxSeqLength + 1];
    std::vector<Matcher::result_t> alnResults;
    alnResults.reserve(maxSetSize);
    std::string backtrace;
    std::string result;

    while (kseq_read(seq) >= 0) {
        if (seq->name.l == 0 || seq->seq.l == 0) {
            Debug(Debug::WARNING) << "Invalid fasta sequence " << setSize << " in entry\n";
            fastaError = true;
            break;
        }

        if (seq->seq.l > maxSeqLength) {
            Debug(Debug::WARNING) << "Member sequence " << setSize << " in entry too long\n";
            fastaError = true;
            break;
        }

        // first sequence is always the query
        if (setSize == 0) {
            centerLengthWithGaps = seq->seq.l;
            backtrace.reserve(centerLengthWithGaps);
            if (maskByFirst == true) {
                for (size_t i = 0; i < centerLengthWithGaps; ++i) {
                    if (seq->seq.s[i] == '-') {
                        maskedColumns[i] = true;
                        maskedCount++;
                    } else {
                        maskedColumns[i] = false;
                    }
                }
            }
        }

        sequence.mapSequence(0, 0, seq->seq.s, seq->seq.l);
        msaSequences[setSize] = msaContent + msaPos;

        for (size_t i = 0; i < centerLengthWithGaps; ++i) {
            if (maskByFirst == true && maskedColumns[i] == true) {
                continue;
            }

            // skip a3m lower letters
            if (msaType == MSA_A3M && islower(seq->seq.s[i])) {
                continue;
            }

            msaContent[msaPos++] = (seq->seq.s[i] == '-') ? (int)MultipleAlignment::GAP : sequence.numSequence[i];
        }

        // fill up the sequence buffer for the SIMD profile calculation
        size_t rowSize = msaPos / (VECSIZE_INT*4);
        rowSize = (rowSize+1) * (VECSIZE_INT*4);
        while(msaPos < rowSize) {
            msaContent[msaPos++] = MultipleAlignment::GAP;
        }

        setSize++;
    }
    
    if (fastaError == true) {
        Debug(Debug::WARNING) << "Invalid msa ! Skipping entry.\n";
        return "";
    }

    if (setSize == 0) {
        Debug(Debug::WARNING) << "Empty msa ! Skipping entry.\n";
        return "";
    }

    if (maskByFirst == false) {

        if (externalMaskedColumns == NULL) {
            PSSMCalculator::computeSequenceWeights(seqWeight, centerLengthWithGaps,
                                                   setSize, const_cast<const char**>(msaSequences));

            // Replace GAP with ENDGAP for all end gaps
            // ENDGAPs are ignored for counting percentage (multi-domain proteins)
            for (unsigned int k = 0; k < setSize; ++k) {
                for (unsigned int i = 0; i < centerLengthWithGaps && msaSequences[k][i] == MultipleAlignment::GAP; ++i)
                    msaSequences[k][i] = MultipleAlignment::ENDGAP;
                for (unsigned int i = centerLengthWithGaps - 1; msaSequences[k][i] == MultipleAlignment::GAP; i--)
                    msaSequences[k][i] = MultipleAlignment::ENDGAP;
            }

            for (unsigned int l = 0; l < centerLengthWithGaps; l++) {
                float res = 0;
                float gap = 0;
                // Add up percentage of gaps
                for (unsigned int k = 0; k < setSize; ++k) {
                    if (msaSequences[k][l] < MultipleAlignment::GAP) {
                        res += seqWeight[k];
                    } else if (msaSequences[k][l] != MultipleAlignment::ENDGAP) {
                        gap += seqWeight[k];
                    } else if (msaSequences[k][l] == MultipleAlignment::ENDGAP) {
                        msaSequences[k][l] = MultipleAlignment::GAP;
                    }
                }

                maskedColumns[l] =  (gap / (res + gap)) > matchRatio;
                maskedCount += maskedColumns[l] ? 1 : 0;
            }

        } else {
            delete[] maskedColumns;
            maskedColumns = externalMaskedColumns;
            for (unsigned int i = 0; i < centerLengthWithGaps; ++i) {
                maskedCount += maskedColumns[i] ? 1 : 0;
            }
        }

        for (unsigned int k = 0; k < setSize; ++k) {
            unsigned int currentCol = 0;
            for (unsigned int l = 0; l < centerLengthWithGaps; ++l) {
                if (maskedColumns[l] == false) {
                    msaSequences[k][currentCol++] = msaSequences[k][l];
                }
            }

            for (unsigned int l = currentCol; l < centerLengthWithGaps; ++l) {
                msaSequences[k][l] = MultipleAlignment::GAP;
            }
        }
    }
    unsigned int centerLength = centerLengthWithGaps - maskedCount;

    MultipleAlignment::MSAResult msaResult(centerLength, centerLength, setSize, msaSequences);
    size_t filteredSetSize = setSize;
    if (filterMsa == 1) {
        filteredSetSize = filter.filter(setSize, centerLength, static_cast<int>(covMSAThr * 100),
                                        qid_vec, qsc,
                                        static_cast<int>(filterMaxSeqId * 100), Ndiff, filterMinEnable,
                                        (const char **) msaSequences, true);
    }

    PSSMCalculator::Profile pssmRes =
            pssmCalculator.computePSSMFromMSA(filteredSetSize, msaResult.centerLength,
                                              (const char **) msaResult.msaSequence, alnResults, wg, scoreBias);
    if (compBiasCorrection == true) {
        SubstitutionMatrix::calcGlobalAaBiasCorrection(&subMat, pssmRes.pssm, pNullBuffer,
                                                       Sequence::PROFILE_AA_SIZE,
                                                       centerLength);
    }
    unsigned char * consensus = new unsigned char[centerLength];
    for (size_t i = 0; i < centerLength; ++i)
        consensus[i] = subMat.aa2num[pssmRes.consensus[i]];
    pssmRes.toBuffer(consensus, centerLength, subMat, result);

    if (externalMaskedColumns == NULL) {
        // Save mask if external mask not given
        result.push_back('\n');
        for (size_t z = 0; z < centerLengthWithGaps; ++z)
            result.push_back(maskedColumns[z] == false ? '0' : '1');
        delete[] maskedColumns;
    }
    delete[] seqWeight;
    delete[] pNullBuffer;
    free(msaSequences);
    free(msaContent);
    // result.push_back('\0');
    
    return result;
}

// Map 0001100 to [ 0 1 2 5 6 ]
// needs to be ungapped->gapped direction
std::vector<int> maskToMapping(std::string mask, size_t length) {
    std::vector<int> mapping;
    mapping.reserve(length); // length of actual sequence
    for (size_t i = 0; i < mask.length(); ++i) {
        if (mask[i] == '0')
            mapping.push_back(i);
    }
    return mapping;
}

enum State {
    SEQ = 0,
    GAP = 1
};

struct Instruction {
    int state;
    int count;
    Instruction(int i_state, int i_count) : state(i_state), count(i_count) {};
    void print() {
        char state_char = (state == SEQ) ? 'S' : 'G';
        std::cout << state_char << " " << count << std::endl;
    }
    char stateChar() { return (state == SEQ) ? 'S' : 'G'; }
};

/**
 * @brief Get merge instructions for two MSAs
 * 
 * @param res  - alignment result
 * @param map1 - ungapped->gapped mapping for msa1
 * @param map2 - ungapped->gapped mapping for msa2
 * @param qBt  - vector to store query merge instructions
 * @param tBt  - vector to store target merge instructions
 */
void getMergeInstructions(
    Matcher::result_t &res,
    std::vector<int> map1,
    std::vector<int> map2,
    std::vector<Instruction> &qBt,
    std::vector<Instruction> &tBt
) {
    qBt.emplace_back(SEQ, 1);  // first match
    tBt.emplace_back(SEQ, 1);
    int new_q, dq;
    int new_t, dt;
    int old_q = map1[res.qStartPos];
    int old_t = map2[res.dbStartPos];
    int q = res.qStartPos + 1;  // indices in non-gappy sequence
    int t = res.dbStartPos + 1;
    
    // Generate instructions for query/target sequences from backtrace
    for (size_t i = 1; i < res.backtrace.length(); ++i) {
        switch (res.backtrace[i]) {
            case 'M': {
                new_q = map1[q];
                new_t = map2[t];
                dq = new_q - old_q;
                dt = new_t - old_t; 
                if (dq == 0) {
                    // No matches in query
                    if (dt > 0)
                        qBt.emplace_back(GAP, dt); 
                    tBt.emplace_back(SEQ, dt);
                } else if (dq == 1) {
                    // One match in query
                    if ((dt - 1) > 0)
                        qBt.emplace_back(GAP, dt - 1);
                    qBt.emplace_back(SEQ, 1);
                    tBt.emplace_back(SEQ, dt);
                } else if (dq >= dt) {
                    // More query matches than target
                    qBt.emplace_back(SEQ, dq);
                    tBt.emplace_back(GAP, dq - dt);
                    tBt.emplace_back(SEQ, dt);
                } else if (dt > dq) {
                    // More target than query
                    qBt.emplace_back(GAP, dt - dq);
                    qBt.emplace_back(SEQ, dq);
                    tBt.emplace_back(SEQ, dt);
                }
                old_q = new_q;
                old_t = new_t;
                ++q;
                ++t;
                break;
            }
            case 'I': {
                ++q;
                break;
            }
            case 'D': {
                ++t;
                break;
            }
        }
    }
}

/**
 * @brief Merges two MSAs
 * 
 * @param msa1 - query MSA
 * @param msa2 - target MSA
 * @param res  - alignment result
 * @param map1 - ungapped->gapped mapping for msa1
 * @param map2 - ungapped->gapped mapping for msa2
 * @param qBt  - query merge instructions
 * @param tBt  - target merge instructions
 * @return std::string - merged MSA
 */
std::string mergeTwoMsa(
    std::string &msa1,
    std::string &msa2,
    Matcher::result_t &res,
    std::vector<int> map1,
    std::vector<int> map2,
    std::vector<Instruction> &qBt,
    std::vector<Instruction> &tBt
) {
    // Calculate pre/end gaps/sequences from backtrace
    size_t qPreSequence = map1[res.qStartPos];
    size_t qPreGaps     = map2[res.dbStartPos];
    size_t qEndSequence = map1[map1.size() - 1] - map1.at(res.qEndPos);
    size_t qEndGaps     = map2[map2.size() - 1] - map2.at(res.dbEndPos);
    size_t tPreSequence = qPreGaps;
    size_t tPreGaps     = qPreSequence;
    size_t tEndSequence = qEndGaps;
    size_t tEndGaps     = qEndSequence;
    
    int q, t;

    // String for merged MSA
    std::string msa; 
    
    // Query msa (msa1) first
    kseq_buffer_t d;
    d.buffer = (char*)msa1.c_str();
    d.length = msa1.size();
    kseq_t *seq = kseq_init(&d);
    while (kseq_read(seq) >= 0) {
        // Header
        msa.push_back('>');
        msa += seq->name.s;
        msa.push_back('\n');
        
        // Pre-alignment: in query, gaps before sequence
        msa.append(qPreGaps, '-');
        msa.append(seq->seq.s, 0, qPreSequence);
        
        // In query, add sequence on M or I, gap on D
        q = qPreSequence;
        for (size_t i = 0; i < qBt.size(); ++i) {
            Instruction ins = qBt[i];
            if (ins.state == SEQ) {
                msa.append(seq->seq.s, q, ins.count);
                q += ins.count;
            } else if (ins.state == GAP) {
                msa.append(ins.count, '-');
            }
        }
        // Post-alignment: in query, sequence before gaps
        msa.append(seq->seq.s, q, qEndSequence);
        msa.append(qEndGaps, '-'); 
        msa.push_back('\n');
    }
    kseq_destroy(seq);
    
    // Target msa (msa2)
    kseq_buffer_t d2;
    d2.buffer = (char*)msa2.c_str();
    d2.length = msa2.size();
    kseq_t *seq2 = kseq_init(&d2);
    while (kseq_read(seq2) >= 0) {
        // Header
        msa.push_back('>');
        msa += seq2->name.s;
        msa.push_back('\n');
        
        // Pre-alignment: in query, gaps before sequence
        msa.append(seq2->seq.s, 0, tPreSequence);
        msa.append(tPreGaps, '-');
        
        // In query, add sequence on M or I, gap on D
        t = tPreSequence;
        for (size_t i = 0; i < tBt.size(); ++i) {
            Instruction ins = tBt[i];
            if (ins.state == SEQ) {
                msa.append(seq2->seq.s, t, ins.count);
                t += ins.count;
            } else if (ins.state == GAP) {
                msa.append(ins.count, '-');
            }
        }
        // Post-alignment: in query, sequence before gaps
        msa.append(tEndGaps, '-');
        msa.append(seq2->seq.s, t, tEndSequence);
        msa.push_back('\n');
    }
    // remove \n
    // msa.erase(msa.length() - 1, 1);
    kseq_destroy(seq2);
    
    return msa;
}

int structuremsa(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    par.compBiasCorrection = 0;

    // Databases
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    DBReader<unsigned int> seqDbrAA(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbrAA.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbr3Di((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbr3Di.open(DBReader<unsigned int>::NOSORT);
   
    IndexReader qdbrH(par.db1, par.threads, IndexReader::HEADERS, touch ? IndexReader::PRELOAD_INDEX : 0);
    
    std::cout << "Got databases" << std::endl;
   
    SubstitutionMatrix subMat_3di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias3di);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    SubstitutionMatrix subMat_aa(blosum.c_str(), 1.4, par.scoreBiasAa);

    std::cout << "Got substitution matrices" << std::endl;
    
    // Initialise MSAs, Sequence objects
    int sequenceCnt = seqDbrAA.getSize();
    std::vector<Sequence*> allSeqs_aa(sequenceCnt);
    std::vector<Sequence*> allSeqs_3di(sequenceCnt);
    std::string msa_aa[sequenceCnt];
    std::string msa_3di[sequenceCnt];
    std::string headers[sequenceCnt];
    std::string mappings[sequenceCnt];
    int idMappings[sequenceCnt];
    std::map<std::string, int> headers_rev;

    int maxSeqLength = 0;
    for (int i = 0; i < sequenceCnt; i++) {
        unsigned int seqKeyAA = seqDbrAA.getDbKey(i);
        unsigned int seqKey3Di = seqDbr3Di.getDbKey(i);
        allSeqs_aa[i] = new Sequence(par.maxSeqLen, seqDbrAA.getDbtype(), (const BaseMatrix *) &subMat_aa, 0, false, par.compBiasCorrection);
        allSeqs_aa[i]->mapSequence(i, seqKeyAA, seqDbrAA.getData(i, 0), seqDbrAA.getSeqLen(i));
        allSeqs_3di[i] = new Sequence(par.maxSeqLen, seqDbr3Di.getDbtype(), (const BaseMatrix *) &subMat_3di, 0, false, par.compBiasCorrection);
        allSeqs_3di[i]->mapSequence(i, seqKey3Di, seqDbr3Di.getData(i, 0), seqDbr3Di.getSeqLen(i));
        maxSeqLength = std::max(maxSeqLength, allSeqs_aa[i]->L);
        msa_aa[i] += ">" + SSTR(i) + "\n";
        msa_aa[i] += seqDbrAA.getData(i, 0);
        msa_3di[i] += ">" +  SSTR(i) + "\n";
        msa_3di[i] += seqDbr3Di.getData(i, 0);
        mappings[i] = std::string(seqDbrAA.getSeqLen(i), '0');
        
        // Map each sequence id to itself for now
        idMappings[i] = i;
        
        // Grab headers, remove \0
        std::string header = qdbrH.sequenceReader->getData(seqKeyAA, 0);
        header = header.substr(0, std::min(header.length() - 1, header.find(' ', 0)));
        headers[i] = header;
        headers_rev[header] = i;
    }
    
    // TODO: dynamically calculate and re-init PSSMCalculator/MsaFilter each iteration
    maxSeqLength = par.maxSeqLen;
    std::cout << "Initialised MSAs, Sequence objects" << std::endl;

    // Substitution matrices needed for query profile
    int8_t *tinySubMatAA  = (int8_t*) mem_align(ALIGN_INT, subMat_aa.alphabetSize * 32);
    int8_t *tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat_3di.alphabetSize * 32);
    for (int i = 0; i < subMat_3di.alphabetSize; i++)
        for (int j = 0; j < subMat_3di.alphabetSize; j++)
            tinySubMat3Di[i * subMat_3di.alphabetSize + j] = subMat_3di.subMatrix[i][j]; // for farrar profile
    for (int i = 0; i < subMat_aa.alphabetSize; i++)
        for (int j = 0; j < subMat_aa.alphabetSize; j++)
            tinySubMatAA[i * subMat_aa.alphabetSize + j] = subMat_aa.subMatrix[i][j];
    std::cout << "Set up tiny substitution matrices" << std::endl;

    bool * alreadyMerged = new bool[sequenceCnt];
    memset(alreadyMerged, 0, sizeof(bool) * sequenceCnt);
    
    std::vector<AlnSimple> hits;
    if (par.guideTree != "") {
        std::cout << "Loading guide tree: " << par.guideTree << "\n"; 
        std::string tree;
        std::string line;
        std::ifstream newick(par.guideTree);
        if (newick.is_open()) {
            while (std::getline(newick, line))
                tree += line;
            newick.close();
        }
        hits = parseNewick(tree, headers_rev);
        if (par.regressive)
            std::reverse(hits.begin(), hits.end());
    } else {
        hits = updateAllScores(
            tinySubMatAA,
            tinySubMat3Di,
            &subMat_aa,
            allSeqs_aa,
            allSeqs_3di,
            alreadyMerged,
            par.maxSeqLen,
            subMat_3di.alphabetSize,
            par.compBiasCorrection,
            par.compBiasCorrectionScale
        );
        sortHitsByScore(hits);
        std::cout << "Performed initial all vs all alignments\n";
        
        hits = mst(hits, sequenceCnt);
        std::cout << "Generated guide tree\n";
    }

    std::vector<int> merges;

    std::cout << "Optimising merge order\n";
    hits = reorderLinkage(hits, merges, sequenceCnt);
    
    int idx = 0;
    for (size_t i = 0; i < merges.size(); i++) {
        std::cout << "Merging " << merges[i] << " sequences\n";
        for (size_t j = 0; j < merges[i]; j++) {
            std::cout << "  " << hits[idx + j].queryId << "\t" << hits[idx + j].targetId << '\n';
        }
        idx += merges[i];
    }

    std::string nw = orderToTree(hits, headers, sequenceCnt);
    std::cout << "Tree: " << nw << ";\n";

    std::cout << "Merging:\n";
    size_t merged = 0;

#pragma omp parallel
{
    // Initialise alignment objects per thread
    StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat_3di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale);
    PSSMCalculator calculator_aa(&subMat_aa, maxSeqLength + 1, sequenceCnt + 1, par.pcmode, par.pcaAa, par.pcbAa, par.gapOpen.values.aminoacid(), par.gapPseudoCount);
    MsaFilter filter_aa(maxSeqLength + 1, sequenceCnt + 1, &subMat_aa, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
    PSSMCalculator calculator_3di(&subMat_3di, maxSeqLength + 1, sequenceCnt + 1, par.pcmode, par.pca3di, par.pcb3di, par.gapOpen.values.aminoacid(), par.gapPseudoCount);
    MsaFilter filter_3di(maxSeqLength + 1, sequenceCnt + 1, &subMat_3di, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid()); 

    int index = 0; // in hit list
    for (size_t i = 0; i < merges.size(); i++) {

#pragma omp for schedule(dynamic, 1)
        for (int j = 0; j < merges[i]; j++) {
            unsigned int mergedId = std::min(hits[index + j].queryId, hits[index + j].targetId);
            unsigned int targetId = std::max(hits[index + j].queryId, hits[index + j].targetId);
            mergedId = idMappings[mergedId];
            targetId = idMappings[targetId];
            bool queryIsProfile = (Parameters::isEqualDbtype(allSeqs_aa[mergedId]->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));
            bool targetIsProfile = (Parameters::isEqualDbtype(allSeqs_aa[targetId]->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));
 
            // Always merge onto sequence with most information
            if (targetIsProfile && !queryIsProfile) {
                std::swap(mergedId, targetId);
            } else if (targetIsProfile && queryIsProfile) {
                float q_neff_sum = 0.0;
                float t_neff_sum = 0.0;
                for (int i = 0; i < allSeqs_3di[mergedId]->L; i++)
                    q_neff_sum += allSeqs_3di[mergedId]->neffM[i];
                for (int i = 0; i < allSeqs_3di[targetId]->L; i++)
                    t_neff_sum += allSeqs_3di[targetId]->neffM[i];
                if (q_neff_sum <= t_neff_sum)
                    std::swap(mergedId, targetId);
            }
 
            assert(mergedId != targetId);

            // Make sure all relevant ids are updated
            for (int k = 0; k < sequenceCnt; k++) {
                if (idMappings[k] == (int)targetId || idMappings[k] == (int)mergedId)
                    idMappings[k] = mergedId;
            }

            structureSmithWaterman.ssw_init(
                allSeqs_aa[mergedId],
                allSeqs_3di[mergedId],
                tinySubMatAA,
                tinySubMat3Di,
                &subMat_aa
            );
            Matcher::result_t res = pairwiseAlignment(
                structureSmithWaterman,
                allSeqs_aa[mergedId]->L,
                allSeqs_aa[targetId],
                allSeqs_3di[targetId],
                par.gapOpen.values.aminoacid(),
                par.gapExtend.values.aminoacid()
            );
            
            // Convert 010101 mask to [ 0, 2, 4 ] index mapping
            std::vector<int> map1 = maskToMapping(mappings[mergedId], res.qLen);
            std::vector<int> map2 = maskToMapping(mappings[targetId], res.dbLen);
            
            // Save new MSAs and remove targetId MSAs
            std::vector<Instruction> qBt;
            std::vector<Instruction> tBt;
            getMergeInstructions(res, map1, map2, qBt, tBt);
            msa_aa[mergedId] = mergeTwoMsa(msa_aa[mergedId], msa_aa[targetId], res, map1, map2, qBt, tBt);
            msa_aa[targetId] = "";
            msa_3di[mergedId] = mergeTwoMsa(msa_3di[mergedId], msa_3di[targetId], res, map1, map2, qBt, tBt);
            msa_3di[targetId] = "";
            assert(msa_aa[mergedId].length() == msa_3di[mergedId].length());

            std::string profile_aa = fastamsa2profile(
                msa_aa[mergedId], calculator_aa, filter_aa, subMat_aa, maxSeqLength,
                sequenceCnt + 1, par.matchRatio, par.filterMsa,
                par.compBiasCorrection,
                par.qid, par.filterMaxSeqId, par.Ndiff, 0,
                par.qsc, par.filterMinEnable, par.wg, NULL, par.scoreBiasAa
            );
            // Mapping is stored at the end of the profile (to \n), so save to mappings[]
            // Iterate backwards until newline to recover the full mask
            std::string mask;
            for (size_t k = profile_aa.length() - 1; profile_aa[k] != '\n'; k--)
                mask.push_back(profile_aa[k]);
            std::reverse(mask.begin(), mask.end());
            mappings[mergedId] = mask;
            
            // Remove mask from the profile itself, -1 for \n
            profile_aa.erase(profile_aa.length() - mappings[mergedId].length() - 1);
            
            // Convert back to bool array to pass to 3di fastamsa2profile
            bool *maskBool = new bool[mask.length()];
            for (size_t k = 0; k < mask.length(); ++k)
                maskBool[k] = (mask[k] == '1') ? true : false;
            
            std::string profile_3di = fastamsa2profile(
                msa_3di[mergedId], calculator_3di, filter_3di, subMat_3di, maxSeqLength,
                sequenceCnt + 1, par.matchRatio, par.filterMsa,
                par.compBiasCorrection,
                par.qid, par.filterMaxSeqId, par.Ndiff, par.covMSAThr,
                par.qsc,
                par.filterMinEnable, par.wg, maskBool, par.scoreBias3di
            );
            delete[] maskBool; 
            assert(profile_aa.length() == profile_3di.length());
            
            if (Parameters::isEqualDbtype(allSeqs_aa[mergedId]->getSeqType(), Parameters::DBTYPE_AMINO_ACIDS)) {
                delete allSeqs_aa[mergedId];
                delete allSeqs_3di[mergedId];
                allSeqs_aa[mergedId] = new Sequence(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_aa, 0, false, par.compBiasCorrection);
                allSeqs_3di[mergedId] = new Sequence(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_3di, 0, false, par.compBiasCorrection);
            }

            allSeqs_aa[mergedId]->mapSequence(mergedId, mergedId, profile_aa.c_str(), profile_aa.length() / Sequence::PROFILE_READIN_SIZE);
            allSeqs_3di[mergedId]->mapSequence(mergedId, mergedId, profile_3di.c_str(), profile_3di.length() / Sequence::PROFILE_READIN_SIZE);
            alreadyMerged[targetId] = true;
        }
        index += merges[i];
        merged += merges[i];
    }
}
    // Cleanup
    seqDbrAA.close();
    seqDbr3Di.close();
    delete[] alreadyMerged;
    delete [] tinySubMatAA;
    delete [] tinySubMat3Di;
    for(size_t i = 0; i < allSeqs_aa.size(); i++){
        delete allSeqs_aa[i];
        delete allSeqs_3di[i];
    }

    // Find the final MSA (only non-empty string left in msa vectors)
    int msaCnt = 0;
    std::string finalMSA;
    for (int i = 0; i < sequenceCnt; ++i) {
        if (msa_aa[i] != "" && msa_3di[i] != "") {
            finalMSA = msa_aa[i];
            ++msaCnt;
            continue;
        }
    }
    // std::cout << finalMSA << '\n';
    assert(msaCnt == 1);
    assert(finalMSA != "");

    // Write final MSA to file with correct headers
    DBWriter resultWriter(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_OMIT_FILE);
    resultWriter.open();
    kseq_buffer_t d;
    d.buffer = (char*)finalMSA.c_str();
    d.length = finalMSA.length();
    kseq_t *seq = kseq_init(&d);
    resultWriter.writeStart(0);
    std::string buffer;
    buffer.reserve(10 * 1024);
    while (kseq_read(seq) >= 0) {
        buffer.append(1, '>');
        buffer.append(headers[std::stoi(seq->name.s)]);
        buffer.append(1, '\n');
        buffer.append(seq->seq.s, seq->seq.l);
        buffer.append(1, '\n');
        resultWriter.writeAdd(buffer.c_str(), buffer.size(), 0);
        buffer.clear();
    }
    resultWriter.writeEnd(0, 0, false, 0);
    resultWriter.close(true);
    FileUtil::remove(par.db2Index.c_str());
   
    return EXIT_SUCCESS;
}