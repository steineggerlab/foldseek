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
#include "KSeqWrapper.h"
#include "LDDT.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "msa2lddt.h"

#include "refinemsa.h"
#include "structuremsa.h"

#ifdef OPENMP
#include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define	EXIT_FAILURE	1
#define	EXIT_SUCCESS	0

struct AlnSimple {
    size_t queryId;
    size_t targetId;
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

// FIXME inconsistent with tree produced by orderToTree
std::vector<AlnSimple> parseNewick(std::string newick, std::map<std::string, int> &headers) {
    // Should know this from number of structures in database
    int total = std::count(newick.begin(), newick.end(), '(');

    // Tokenize tree on ; | ( ) , :
    // Use {-1, 0} on sregex_token_iterator to get matches AND inbetween (i.e. names)
    std::regex pattern("\\s*(;|\\(|\\)|,|:)\\s*");
    auto words_begin = std::sregex_token_iterator(newick.begin(), newick.end(), pattern, {-1, 0});
    auto words_end = std::sregex_token_iterator();
    
    // Initialise TNode array (2n+1 but a +1 for the root)
    // TNode nodes[total * 2 + 2];
    std::vector<TNode *> nodes(total * 2 + 2);
    std::vector<TNode *> parents;
    TNode *tree;
    TNode *subtree;
    
    // Initialise the root node (the +1)
    nodes[0] = new TNode();
    nodes[0]->id = 0;
    tree = nodes[0];
    
    int count = 1;
    std::string prevToken;
    
    for (std::sregex_token_iterator i = words_begin; i != words_end; ++i) {
        std::string match_str = *i;

        if (match_str == "(") {
            // add new node, set it as new subtree
            nodes[count] = new TNode();
            nodes[count]->id = count;
            nodes[count]->dbId = -1;
            nodes[count]->length = 1.0;
            subtree = nodes[count];
            count++;
            
            // add it as child to current tree
            tree->children.push_back(subtree);
            subtree->parent = tree;
            
            // add the tree as parent, set subtree as tree
            parents.push_back(tree);
            tree = subtree;
        } else if (match_str == ",") {
            nodes[count] = new TNode();
            nodes[count]->id = count;
            nodes[count]->dbId = -1;
            nodes[count]->length = 1.0;
            subtree = nodes[count];
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
    
    // Cleanup 
    for(size_t i = 0; i < nodes.size(); i++)
        delete nodes[i];

    return hits;
}

int rescoreBacktrace(
    int qpos, int dbpos,
    Sequence *qSeqAA, Sequence *tSeqAA,
    Sequence *qSeq3Di, Sequence *tSeq3Di,
    int gapOpen, int gapExtend, std::string backtrace,
    SubstitutionMatrix *mat_aa,
    SubstitutionMatrix *mat_3di
    // short ** mat_aa, short ** mat_3di
) {
    unsigned int gapIcount, gapDcount;
    int rescore, gapScore;
    
    rescore = 0;

    for(size_t j = 0; j < backtrace.length(); j++){

        //if match add 3Di and AA score as rescore value
        int qAALetter = qSeqAA->numSequence[qpos];
        int dbAALetter =  tSeqAA->numSequence[dbpos];
        int q3DiLetter = qSeq3Di->numSequence[qpos];
        int db3DiLetter =  tSeq3Di->numSequence[dbpos];

        switch (backtrace[j]) {
            case 'M':
                rescore += mat_aa->subMatrix[qAALetter][dbAALetter] + mat_3di->subMatrix[q3DiLetter][db3DiLetter];
                gapIcount = 0; gapDcount = 0;
                qpos++;
                dbpos++;
                break;
            case 'D':
                gapIcount = 0;
                gapScore = (gapDcount == 0) ? -gapOpen : -gapExtend;
                rescore += gapScore;
                gapDcount++;
                dbpos++;
                break;
            case 'I':
                gapDcount = 0;
                gapScore = (gapIcount == 0) ? -gapOpen : -gapExtend;
                rescore += gapScore;
                gapIcount++;
                qpos++;
                break;
        }
    }
    return rescore;
}


Matcher::result_t pairwiseAlignment(
    StructureSmithWaterman & aligner,
    unsigned int querySeqLen,
    Sequence *query_aa,
    Sequence *query_3di,
    Sequence *target_aa,
    Sequence *target_3di,
    int gapOpen, int gapExtend,
    SubstitutionMatrix *mat_aa,
    SubstitutionMatrix *mat_3di,
    std::vector<std::vector<std::vector<int> > > &neighbours,
    std::vector<int> &qMap,
    std::vector<int> &tMap
) {
    std::string backtrace;

    bool targetIsProfile = (Parameters::isEqualDbtype(target_aa->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));
    bool queryIsProfile = (Parameters::isEqualDbtype(query_aa->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));

    unsigned char * query_aa_seq = query_aa->numSequence;
    unsigned char * query_3di_seq = query_3di->numSequence;
    unsigned char * target_aa_seq = target_aa->numSequence;
    unsigned char * target_3di_seq = target_3di->numSequence;
    if (queryIsProfile) {
        query_aa_seq = query_aa->numConsensusSequence;
        query_3di_seq = query_3di->numConsensusSequence;
    }
    if (targetIsProfile) {
        target_aa_seq = target_aa->numConsensusSequence;
        target_3di_seq = target_3di->numConsensusSequence;
    }

    // TODO composition bias
    float *composition_bias_aa  = new float[query_aa->L];
    float *composition_bias_ss  = new float[query_aa->L];
    float *tmp_composition_bias = new float[query_aa->L];
    if (true) {
        SubstitutionMatrix::calcLocalAaBiasCorrection(mat_aa, query_aa->numSequence, query_aa->L, tmp_composition_bias, 1.0);
        for (int i =0; i < query_aa->L; i++) {
            composition_bias_aa[i] = (int8_t) (tmp_composition_bias[i] < 0.0) ? tmp_composition_bias[i] - 0.5 : tmp_composition_bias[i] + 0.5;
        }
        SubstitutionMatrix::calcLocalAaBiasCorrection(mat_3di, query_3di->numSequence, query_3di->L, tmp_composition_bias, 1.0);
        for (int i =0; i < query_aa->L; i++) {
            composition_bias_ss[i] = (int8_t) (tmp_composition_bias[i] < 0.0) ? tmp_composition_bias[i] - 0.5 : tmp_composition_bias[i] + 0.5;
        }
    } else {
        memset(composition_bias_aa, 0, query_aa->L * sizeof(int8_t));
        memset(composition_bias_ss, 0, query_aa->L * sizeof(int8_t));
    }
        
    // StructureSmithWaterman::s_align align = aligner.alignScoreEndPos<StructureSmithWaterman::PROFILE>(
    //     target_aa_seq,
    //     target_3di_seq,
    //     target_aa->L,
    //     gapOpen,
    //     gapExtend,
    //     querySeqLen / 2
    // );
    // align = aligner.alignStartPosBacktrace<StructureSmithWaterman::PROFILE>(
    //     target_aa_seq,
    //     target_3di_seq,
    //     target_aa->L,
    //     gapOpen,
    //     gapExtend,
    //     3,
    //     backtrace,
    //     align,
    //     0,
    //     0.0,
    //     querySeqLen / 2
    // );
    // unsigned int alnLength = Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1, align.dbStartPos1, align.dbEndPos1);
    // alnLength = backtrace.size();
    // float seqId = Util::computeSeqId(Parameters::SEQ_ID_ALN_LEN, align.identicalAACnt, querySeqLen, target_aa->L, alnLength); 
    // Matcher::result_t sw_align(
    //     target_aa->getDbKey(),
    //     align.score1,
    //     align.qCov,
    //     align.tCov,
    //     seqId,
    //     align.evalue,
    //     alnLength,
    //     align.qStartPos1,
    //     align.qEndPos1,
    //     querySeqLen,
    //     align.dbStartPos1,
    //     align.dbEndPos1,
    //     target_aa->L,
    //     backtrace
    // );

    short **query_profile_scores_aa = new short * [aligner.get_profile()->alphabetSize];
    short **query_profile_scores_3di = new short * [aligner.get_profile()->alphabetSize];
    for (int32_t j = 0; j < aligner.get_profile()->alphabetSize; j++) {
        query_profile_scores_aa[j] = new short [querySeqLen];
        query_profile_scores_3di[j] = new short [querySeqLen];
    }
    if (queryIsProfile) {
        for (unsigned int i = 0; i < querySeqLen; i++) {
            for (int32_t j = 0; j < aligner.get_profile()->alphabetSize; j++) {
                query_profile_scores_aa[j][i]  = query_aa->profile_for_alignment[j * querySeqLen + i];
                query_profile_scores_3di[j][i] = query_3di->profile_for_alignment[j * querySeqLen + i];
            }
        }
    } else {
        for (unsigned int i = 0; i < querySeqLen; i++) {
            for (int32_t j = 0; j < aligner.get_profile()->alphabetSize; j++) {
                query_profile_scores_aa[j][i]  = mat_aa->subMatrix[j][query_aa_seq[i]]   + composition_bias_aa[i];
                query_profile_scores_3di[j][i] = mat_3di->subMatrix[j][query_3di_seq[i]] + composition_bias_ss[i];
            }
        }
    }
   
    short **target_profile_scores_aa = new short * [aligner.get_profile()->alphabetSize];
    short **target_profile_scores_3di = new short * [aligner.get_profile()->alphabetSize];
    for (int32_t j = 0; j < aligner.get_profile()->alphabetSize; j++) {
        target_profile_scores_aa[j]  = new short [target_aa->L];
        target_profile_scores_3di[j] = new short [target_aa->L];
    }
    if (targetIsProfile) {
        for (int i = 0; i < target_aa->L; i++) {
            for (int32_t j = 0; j < aligner.get_profile()->alphabetSize; j++) {
                target_profile_scores_aa[j][i]  = target_aa->profile_for_alignment[j * target_aa->L + i];
                target_profile_scores_3di[j][i] = target_3di->profile_for_alignment[j * target_aa->L + i];
            }
        }
    } else {
        for (int i = 0; i < target_aa->L; i++) {
            for (int32_t j = 0; j < aligner.get_profile()->alphabetSize; j++) {
                target_profile_scores_aa[j][i]  = mat_aa->subMatrix[j][target_aa_seq[i]];
                target_profile_scores_3di[j][i] = mat_3di->subMatrix[j][target_3di_seq[i]];
            }
        }
    }

    delete[] composition_bias_aa;
    delete[] composition_bias_ss;
    delete[] tmp_composition_bias;
    
    // for (unsigned int i = 0; i < aligner.get_profile()->alphabetSize; i++) {
    //     for (unsigned int j = 0; j < aligner.get_profile()->alphabetSize; j++) {
    //         std::cout << std::fixed << std::setprecision(3) << query_subMat_psp_3di[i][j] << '\t';
    //     }
    //     std::cout << '\n';
    // }
    // std::cout << '\n';

    Matcher::result_t gAlign = aligner.simpleGotoh(
        target_aa_seq,
        target_3di_seq,
        query_profile_scores_aa,
        query_profile_scores_3di,
        target_profile_scores_aa,
        target_profile_scores_3di,
        0,
        query_aa->L,
        0,
        target_aa->L,
        gapOpen,
        gapExtend,
        targetIsProfile,
        neighbours,
        query_aa->getId(),
        target_aa->getId(),
        qMap,
        tMap
    );
    
    // int sw_rescore = rescoreBacktrace(
    //     sw_align.qStartPos,
    //     sw_align.dbStartPos,
    //     query_aa,
    //     target_aa,
    //     query_3di,
    //     target_3di,
    //     gapOpen, gapExtend,
    //     sw_align.backtrace,
    //     mat_aa,
    //     mat_3di
    // ); 

    // int gt_rescore = rescoreBacktrace(
    //     gAlign.qStartPos,
    //     gAlign.dbStartPos,
    //     query_aa,
    //     target_aa,
    //     query_3di,
    //     target_3di,
    //     gapOpen, gapExtend,
    //     gAlign.backtrace,
    //     mat_aa,
    //     mat_3di
    // ); 

    // std::cout << "\nTarget is profile " << targetIsProfile << ", query " << queryIsProfile << '\n';
    // std::cout << "sw rescore: " << sw_rescore << '\n';
    // std::cout << "gotoh rescore: " << gt_rescore << "\n";
    // std::cout << sw_align.backtrace << '\n' << gAlign.backtrace << '\n';
    // std::cout << " qStartPos: " << sw_align.qStartPos << ", " << gAlign.qStartPos << '\n';
    // std::cout << "   qEndPos: " << sw_align.qEndPos << ", " << gAlign.qEndPos << '\n';
    // std::cout << "dbStartPos: " << sw_align.dbStartPos << ", " << gAlign.dbStartPos << '\n';
    // std::cout << "  dbEndPos: " << sw_align.dbEndPos << ", " << gAlign.dbEndPos << '\n';
    // std::cout << "     Score: " << sw_align.score << ", " << gAlign.score << '\n';

    // assert(sw_align.qStartPos == gAlign.qStartPos);
    // assert(sw_align.dbStartPos == gAlign.dbStartPos);
    // assert(sw_align.qEndPos == gAlign.qEndPos);
    // assert(sw_align.dbEndPos == gAlign.dbEndPos);
    // assert(sw_align.score == gAlign.score);
    // assert(backtrace == gAlign.backtrace);
    // assert(sw_rescore == gt_rescore);

    // if (targetIsProfile && queryIsProfile) {
    //     exit(1);
    // }
    
    for (int32_t i = 0; i < aligner.get_profile()->alphabetSize; i++) {
        delete[] query_profile_scores_aa[i];
        delete[] query_profile_scores_3di[i];
        delete[] target_profile_scores_aa[i];
        delete[] target_profile_scores_3di[i];
    }
    
    delete[] query_profile_scores_aa;
    delete[] query_profile_scores_3di;
    delete[] target_profile_scores_aa;
    delete[] target_profile_scores_3di;
    
    return gAlign;
    // return sw_align;
}

void sortHitsByScore(std::vector<AlnSimple> &hits) {
    SORT_PARALLEL(hits.begin(), hits.end(), [](const AlnSimple & a, const AlnSimple & b) {
        // sort by score then qId then tId
        if (a.score == b.score) {
            if (a.queryId == b.queryId) {
                return a.targetId < b.targetId;
            }
            return a.queryId < b.queryId;
        }
        return a.score > b.score;
    });
}

std::vector<AlnSimple> removeMergedHits(std::vector<AlnSimple> & hits, size_t mergedId, size_t targetId) {
    std::vector<AlnSimple> newHits;
    for (size_t i = 0; i < hits.size(); i++) {
        if (hits[i].queryId != mergedId && hits[i].targetId != mergedId
            && hits[i].queryId != targetId && hits[i].targetId != targetId) {
            newHits.push_back(hits[i]);
        }
    }
    return newHits;
}


inline size_t get1dIndex(size_t i, size_t j, size_t N) {
    return j + i * (2 * N - i - 1) / 2 - i - 1;
}

std::vector<AlnSimple> updateAllScores(
    int8_t * tinySubMatAA,
    int8_t * tinySubMat3Di,
    SubstitutionMatrix * subMat_aa,
    SubstitutionMatrix * subMat_3di,
    std::vector<Sequence*> & allSeqs_aa,
    std::vector<Sequence*> & allSeqs_3di,
    bool * alreadyMerged,
    int maxSeqLen,
    int alphabetSize,
    int compBiasCorrection,
    int compBiasCorrectionScale
) {
    std::vector<AlnSimple> newHits;
#pragma omp parallel
{
    StructureSmithWaterman structureSmithWaterman(
        maxSeqLen,
        alphabetSize,
        compBiasCorrection,
        compBiasCorrectionScale,
        subMat_aa,
        subMat_3di
    );
    std::vector<AlnSimple> threadHits;

#pragma omp for schedule(dynamic, 10)
    for (unsigned int i = 0; i < allSeqs_aa.size(); i++) {
        if (alreadyMerged[i])
            continue;
        structureSmithWaterman.ssw_init(
            allSeqs_aa[i],
            allSeqs_3di[i],
            tinySubMatAA,
            tinySubMat3Di,
            subMat_aa
        );
        for (size_t j = i + 1; j < allSeqs_aa.size(); j++) {
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
            aln.queryId = allSeqs_aa[i]->getId();
            aln.targetId = allSeqs_aa[j]->getId();
            aln.score = structureSmithWaterman.ungapped_alignment(target_aa_seq, target_3di_seq, allSeqs_aa[j]->L);
            threadHits.push_back(aln); 
        }
    }

#pragma omp critical
    {
        newHits.insert(newHits.end(), threadHits.begin(), threadHits.end());
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
std::string orderToTree(std::vector<AlnSimple> hits, std::vector<std::string> &headers, int n) {
    std::string nwk = "";
    int totalNodes = n * 2 + 2;

    std::vector<TNode *> nodes(totalNodes);
    for (int i = 0; i < totalNodes; i++) {
        nodes[i] = new TNode();
        nodes[i]->id = i;
        nodes[i]->parent = nodes[i];
        if (i < n)
            nodes[i]->name = headers[i];
    }
    TNode *root = nodes[0];

    int newId = n + 1;
    for (AlnSimple aln : hits) {
        TNode *u = findRoot(nodes[aln.queryId]);
        TNode *v = findRoot(nodes[aln.targetId]);
        TNode *newNode = nodes[newId];
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
    std::string newick = getNewick(root);
    
    // Cleanup 
    for(size_t i = 0; i < nodes.size(); i++)
        delete nodes[i];

    return newick;
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
std::vector<AlnSimple> reorderLinkage(std::vector<AlnSimple> linkage, std::vector<size_t> &merges, int n) {
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
                                              (const char **) msaResult.msaSequence,
#ifdef GAP_POS_SCORING
                                              alnResults,
#endif
                                              wg);
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
std::vector<int> maskToMapping(std::string mask) {
    std::vector<int> mapping;
    for (size_t i = 0; i < mask.length(); ++i) {
        if (mask[i] == '0')
            mapping.push_back(i);
    }
    return mapping;
}

std::vector<AlnSimple> parseAndScoreExternalHits(
    DBReader<unsigned int> * cluDbr,
    int8_t * tinySubMatAA,
    int8_t * tinySubMat3Di,
    SubstitutionMatrix * subMat_aa,
    SubstitutionMatrix * subMat_3di,
    std::vector<Sequence*> & allSeqs_aa,
    std::vector<Sequence*> & allSeqs_3di,
    int maxSeqLen,
    int alphabetSize,
    int compBiasCorrection,
    int compBiasCorrectionScale
) {
    // open an alignment DBReader
    std::vector<AlnSimple> allAlnResults;

#pragma omp parallel
    {
        StructureSmithWaterman structureSmithWaterman(
            maxSeqLen,
            alphabetSize,
            compBiasCorrection,
            compBiasCorrectionScale,
            subMat_aa,
            subMat_3di
        );
        std::vector<AlnSimple> threadAlnResults;
        char buffer[255 + 1];

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < cluDbr->getSize(); ++i) {
            char *data = cluDbr->getData(i, 0);
            unsigned int queryKey = cluDbr->getDbKey(i);
            structureSmithWaterman.ssw_init(
                allSeqs_aa[queryKey],
                allSeqs_3di[queryKey],
                tinySubMatAA,
                tinySubMat3Di,
                subMat_aa
            );
            while (*data != '\0') {
                Util::parseKey(data, buffer);
                const unsigned int dbKey = (unsigned int) strtoul(buffer, NULL, 10);
                if (queryKey == dbKey) {
                    data = Util::skipLine(data);
                    continue;
                }
                AlnSimple aln;
                aln.queryId = queryKey;
                aln.targetId = dbKey;
                aln.score = structureSmithWaterman.ungapped_alignment(allSeqs_aa[dbKey]->numSequence,
                                                                      allSeqs_3di[dbKey]->numSequence,
                                                                      allSeqs_aa[dbKey]->L);
                threadAlnResults.push_back(aln);
                data = Util::skipLine(data);
            }
        }
#pragma omp critical
        {
            allAlnResults.insert(allAlnResults.end(), threadAlnResults.begin(), threadAlnResults.end());
        }
    }
    return allAlnResults;
}

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
        for (Instruction ins : qBt) {
            if (ins.state == SEQ) {
                msa.append(seq->seq.s, q, ins.count);
                q += ins.count;
            } else if (ins.state == GAP) {
                msa.append(ins.count, '-');
            }
        }
        // Post-alignment: in query, sequence before gaps
        // if (qEndSequence > 0)
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
        if (tEndSequence > 0)
            msa.append(seq2->seq.s, t, tEndSequence);
        msa.push_back('\n');
    }
    // remove \n
    // msa.erase(msa.length() - 1, 1);
    kseq_destroy(seq2);
    
    return msa;
}

void testSeqLens(std::string msa, std::map<std::string, int> &lengths) {
    KSeqWrapper* kseq = new KSeqBuffer(msa.c_str(), msa.length());
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq->entry;
        int count = 0;
        for (size_t i = 0; i < entry.sequence.l; i++) {
            if (entry.sequence.s[i] != '-')
                count++;
        }
        assert(lengths[entry.name.s] == count);
    }
}

void updateBins(
    size_t mergedId,
    size_t targetId,
    Matcher::result_t res,
    std::vector<int> &map1,
    std::vector<int> &map2,
    std::vector<Instruction> &qBt,
    std::vector<Instruction> &tBt,
    std::vector<std::vector<std::vector<int> > > &neighbours
) {
    std::vector<bool> matches;
    size_t qPreSequence = map1[res.qStartPos];
    size_t qPreGaps     = map2[res.dbStartPos];
    size_t qEndGaps     = map2[map2.size() - 1] - map2[res.dbEndPos];
    for (Instruction ins : qBt)
        matches.insert(matches.end(), ins.count, (ins.state == SEQ) ? true : false);
    int m = 0;
    int q = qPreSequence;
    int t = qPreGaps;
    for (Instruction ins : tBt) {
        for (int i = 0; i < ins.count; i++) {
            if (ins.state == SEQ) {
                if (matches[m])
                    for (size_t j = 0; j < neighbours[mergedId][0].size(); j++)
                        neighbours[mergedId][q][j] = std::sqrt(neighbours[mergedId][q][j] * neighbours[targetId][t][j]);
                else
                    neighbours[mergedId].insert(neighbours[mergedId].begin() + q, neighbours[targetId][t]);
                t++;
            }
            m++;
            q++;
        }
    }
    neighbours[mergedId].insert(neighbours[mergedId].begin(), neighbours[targetId].begin(), neighbours[targetId].begin() + qPreGaps);
    neighbours[mergedId].insert(neighbours[mergedId].end(), neighbours[targetId].end() - qEndGaps, neighbours[targetId].end());
}

Matcher::result_t pairwiseTMAlign(
    int mergedId,
    int targetId,
    DBReader<unsigned int> &seqDbrAA,
    DBReader<unsigned int> &seqDbrCA
) {
    int qLen = seqDbrAA.getSeqLen(mergedId);
    int tLen = seqDbrAA.getSeqLen(targetId);
    
    unsigned int qKey = seqDbrAA.getDbKey(mergedId);
    size_t qCaId = seqDbrCA.getId(qKey);

    unsigned int tKey = seqDbrAA.getDbKey(targetId);
    size_t tCaId = seqDbrCA.getId(tKey);
    
    Coordinate16 qcoords;

    char *qcadata = seqDbrCA.getData(qCaId, 0);
    size_t qCaLength = seqDbrCA.getEntryLen(qCaId);
    float *qCaData = qcoords.read(qcadata, qLen, qCaLength);
    char *merged_aa_seq = seqDbrAA.getData(qCaId, 0);
    
    Coordinate16 tcoords;
    char *tcadata = seqDbrCA.getData(tCaId, 0);
    size_t tCaLength = seqDbrCA.getEntryLen(tCaId);
    float *tCaData = tcoords.read(tcadata, tLen, tCaLength);
    char *target_aa_seq = seqDbrAA.getData(tCaId, 0);

    float TMscore = 0.0;
    TMaligner tmaln(std::max(qLen, tLen)+VECSIZE_FLOAT, 1, 0);
    tmaln.initQuery(qCaData, &qCaData[qLen], &qCaData[qLen * 2], merged_aa_seq, qLen);
    Matcher::result_t res = tmaln.align(targetId, tCaData, &tCaData[tLen], &tCaData[tLen * 2], target_aa_seq, tLen, TMscore);
    res.backtrace = Matcher::uncompressAlignment(res.backtrace);
    res.score /= 100;

    return res;
}


int structuremsa(int argc, const char **argv, const Command& command, bool preCluster) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    // Databases
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> seqDbrAA(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_LOOKUP_REV);
    seqDbrAA.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbr3Di((par.db1+"_ss").c_str(), (par.db1+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbr3Di.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> seqDbrCA((par.db1+"_ca").c_str(), (par.db1+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbrCA.open(DBReader<unsigned int>::NOSORT);
   
    IndexReader qdbrH(par.db1, par.threads, IndexReader::HEADERS, touch ? IndexReader::PRELOAD_INDEX : 0);
    
    std::cout << "Got databases" << std::endl;
   
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

    std::cout << "Got substitution matrices" << std::endl;

    // Initialise MSAs, Sequence objects
    size_t sequenceCnt = seqDbrAA.getSize();
    std::vector<Sequence*> allSeqs_aa(sequenceCnt);
    std::vector<Sequence*> allSeqs_3di(sequenceCnt);
    std::vector<std::string> msa_aa(sequenceCnt);
    std::vector<std::string> msa_3di(sequenceCnt);
    std::vector<std::string> headers(sequenceCnt);
    std::vector<std::string> mappings(sequenceCnt);
    std::vector<size_t> idMappings(sequenceCnt);
    std::map<std::string, int> headers_rev;

    std::map<std::string, int> seqLens;

    // 6x6 neighbour matrices
    std::vector<std::vector<std::vector<int> > > neighbours(sequenceCnt);

    int maxSeqLength = par.maxSeqLen;
    for (size_t i = 0; i < sequenceCnt; i++) {
        size_t seqKeyAA = seqDbrAA.getDbKey(i);
        size_t seqKey3Di = seqDbr3Di.getDbKey(i);

        // Grab headers, remove \0
        std::string header = qdbrH.sequenceReader->getData(seqKeyAA, 0);
        header = header.substr(0, std::min(header.length() - 1, header.find(' ', 0)));
        headers[i] = header;
        headers_rev[header] = i;

        // Create Sequences
        allSeqs_aa[i] = new Sequence(par.maxSeqLen, seqDbrAA.getDbtype(), (const BaseMatrix *) &subMat_aa, 0, false, par.compBiasCorrection);
        allSeqs_aa[i]->mapSequence(i, seqKeyAA, seqDbrAA.getData(i, 0), seqDbrAA.getSeqLen(i));
        allSeqs_3di[i] = new Sequence(par.maxSeqLen, seqDbr3Di.getDbtype(), (const BaseMatrix *) &subMat_3di, 0, false, par.compBiasCorrection);
        allSeqs_3di[i]->mapSequence(i, seqKey3Di, seqDbr3Di.getData(i, 0), seqDbr3Di.getSeqLen(i));
        
        maxSeqLength = std::max(maxSeqLength, allSeqs_aa[i]->L);
        msa_aa[i] += ">" + header + "\n";
        msa_aa[i] += seqDbrAA.getData(i, 0);
        msa_3di[i] += ">" +  header + "\n";
        msa_3di[i] += seqDbr3Di.getData(i, 0);
        mappings[i] = std::string(seqDbrAA.getSeqLen(i), '0');

        // Map each sequence id to itself for now
        idMappings[i] = i;
        
        seqLens[header] = allSeqs_3di[i]->L;
        
        Coordinate16 coords;
        char *cadata = seqDbrCA.getData(i, 0);
        size_t caLength = seqDbrCA.getEntryLen(i);
        float *caData = coords.read(cadata, allSeqs_aa[i]->L, caLength);
        neighbours[i].resize(allSeqs_aa[i]->L);
        for (size_t j = 0; j < allSeqs_aa[i]->L; j++) {
            neighbours[i][j].resize(8);
            for (size_t k = 0; k < allSeqs_aa[i]->L; k++) {
                if (j == k) {
                    continue;
                } else {
                    // 5 bins, each 4 angstrom
                    // could store as single integer, each 4 bits
                    double dX = pow(caData[k] - caData[j], 2);
                    double dY = pow(caData[k + allSeqs_aa[i]->L] - caData[j + allSeqs_aa[i]->L], 2);
                    double dZ = pow(caData[k + allSeqs_aa[i]->L * 2] - caData[j + allSeqs_aa[i]->L * 2], 2);
                    double d  = sqrt(dX + dY + dZ);
                    if (d < 5) {
                        neighbours[i][j][0] += d;
                    } else if (d < 7) {
                        neighbours[i][j][1] += d;
                    } else if (d < 9) {
                        neighbours[i][j][2] += d;
                    } else if (d < 10) {
                        neighbours[i][j][3] += d;
                    } else if (d < 11) {
                        neighbours[i][j][4] += d;
                    } else if (d < 12) {
                        neighbours[i][j][5] += d;
                    } else if (d < 13) {
                        neighbours[i][j][6] += d;
                    } else if (d < 15) {
                        neighbours[i][j][7] += d;
                    }
                }
            }
        }
        // for (size_t j = 0; j < allSeqs_aa[i]->L; j++) {
        //     std::cout << headers[i] << '\t'
        //         << neighbours[i][j][0] << '\t'
        //         << neighbours[i][j][1] << '\t'
        //         << neighbours[i][j][2] << '\t'
        //         << neighbours[i][j][3] << '\t'
        //         << neighbours[i][j][4]
        //         << '\n';
        // }
    }
    
    // TODO: dynamically calculate and re-init PSSMCalculator/MsaFilter each iteration
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
   
    DBReader<unsigned int> * cluDbr = NULL;

    if (preCluster) {
        // consider everything merged and unmerge the ones that are not
        memset(alreadyMerged, 1, sizeof(bool) * sequenceCnt);
        cluDbr = new DBReader<unsigned int>(
            par.db2.c_str(),
            par.db2Index.c_str(),
            par.threads,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA
        );
        cluDbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
        // mark all sequences that are already clustered as merged
        for(size_t i = 0; i < cluDbr->getSize(); i++){
            unsigned int dbKey = cluDbr->getDbKey(i);
            alreadyMerged[dbKey] = 0;
        }
    } else {
        memset(alreadyMerged, 0, sizeof(bool) * sequenceCnt);
    }       

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
            &subMat_3di,
            allSeqs_aa,
            allSeqs_3di,
            alreadyMerged,
            par.maxSeqLen,
            subMat_3di.alphabetSize,
            par.compBiasCorrection,
            par.compBiasCorrectionScale
        );
        if (cluDbr != NULL) {
            // add external hits to the list
            std::vector<AlnSimple> externalHits = parseAndScoreExternalHits(
                cluDbr,
                tinySubMatAA,
                tinySubMat3Di,
                &subMat_aa,
                &subMat_3di,
                allSeqs_aa,
                allSeqs_3di,
                par.maxSeqLen,
                subMat_3di.alphabetSize,
                par.compBiasCorrection,
                par.compBiasCorrectionScale
            );
            // maybe a bit dangerous because memory of hits might be doubled
            for (size_t i = 0; i < externalHits.size(); i++)
                hits.push_back(externalHits[i]);
        }
        sortHitsByScore(hits);
        std::cout << "Performed initial all vs all alignments\n";
        
        hits = mst(hits, sequenceCnt);
        std::cout << "Generated guide tree\n";
    }

    std::cout << "Optimising merge order\n";
    std::vector<size_t> merges;
    hits = reorderLinkage(hits, merges, sequenceCnt);

    int idx = 0;
    for (size_t i = 0; i < merges.size(); i++) {
        std::cout << "Merging " << merges[i] << " sequences\n";
        for (size_t j = 0; j < merges[i]; j++) {
            std::cout << "  " << headers[hits[idx + j].queryId] << "\t" << headers[hits[idx + j].targetId] << '\n';
        }
        idx += merges[i];
    }
    
    int finalMSAId = 0;
    std::string finalMSA_aa;
    std::string finalMSA_3di;

    std::string nw = orderToTree(hits, headers, sequenceCnt);
    std::cout << "Tree: " << nw << ";\n";

    std::cout << "Merging:\n";

    int numMerges = 0;

#pragma omp parallel
{
    // Initialise alignment objects per thread
    StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat_3di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMat_aa, &subMat_3di);
    PSSMCalculator calculator_aa(&subMat_aa, maxSeqLength + 1, sequenceCnt + 1, par.pcmode, par.pcaAa, par.pcbAa
#ifdef GAP_POS_SCORING
    , par.gapOpen.values.aminoacid(), par.gapPseudoCount
#endif
    );
    MsaFilter filter_aa(maxSeqLength + 1, sequenceCnt + 1, &subMat_aa, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
    PSSMCalculator calculator_3di(&subMat_3di, maxSeqLength + 1, sequenceCnt + 1, par.pcmode, par.pca3di, par.pcb3di
#ifdef GAP_POS_SCORING
    , par.gapOpen.values.aminoacid(), par.gapPseudoCount
#endif
    );
    MsaFilter filter_3di(maxSeqLength + 1, sequenceCnt + 1, &subMat_3di, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid()); 

    int index = 0; // in hit list
    for (size_t i = 0; i < merges.size(); i++) {

#pragma omp for schedule(dynamic, 1)
        for (size_t j = 0; j < merges[i]; j++) {
            size_t mergedId = std::min(hits[index + j].queryId, hits[index + j].targetId);
            size_t targetId = std::max(hits[index + j].queryId, hits[index + j].targetId);
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

            queryIsProfile = (Parameters::isEqualDbtype(allSeqs_aa[mergedId]->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));
            targetIsProfile = (Parameters::isEqualDbtype(allSeqs_aa[targetId]->getSeqType(), Parameters::DBTYPE_HMM_PROFILE));

            assert(mergedId != targetId);
            finalMSAId = mergedId;

            // Make sure all relevant ids are updated
            for (size_t k = 0; k < sequenceCnt; k++) {
                if (idMappings[k] == targetId || idMappings[k] == mergedId)
                    idMappings[k] = mergedId;
            }

            // Convert 010101 mask to [ 0, 2, 4 ] index mapping
            std::vector<int> map1 = maskToMapping(mappings[mergedId]);
            std::vector<int> map2 = maskToMapping(mappings[targetId]);
            
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
                allSeqs_aa[mergedId],
                allSeqs_3di[mergedId],
                allSeqs_aa[targetId],
                allSeqs_3di[targetId],
                par.gapOpen.values.aminoacid(),
                par.gapExtend.values.aminoacid(),
                &subMat_aa,
                &subMat_3di,
                neighbours,
                map1,
                map2
            );
            std::vector<Instruction> qBt;
            std::vector<Instruction> tBt;
            getMergeInstructions(res, map1, map2, qBt, tBt);
            std::string msa3di_aa  = mergeTwoMsa(msa_aa[mergedId], msa_aa[targetId], res, map1, map2, qBt, tBt);
            std::string msa3di_3di = mergeTwoMsa(msa_3di[mergedId], msa_3di[targetId], res, map1, map2, qBt, tBt);

            // If neither are profiles, do TM-align as well and take the best alignment
            if (!queryIsProfile && !targetIsProfile) {
                Matcher::result_t tmRes = pairwiseTMAlign(mergedId, targetId, seqDbrAA, seqDbrCA);
                std::vector<Instruction> qBtTM;
                std::vector<Instruction> tBtTM;
                getMergeInstructions(tmRes, map1, map2, qBtTM, tBtTM);
                std::string msaTM_aa  = mergeTwoMsa(msa_aa[mergedId],  msa_aa[targetId],  tmRes, map1, map2, qBtTM, tBtTM);
                std::string msaTM_3di = mergeTwoMsa(msa_3di[mergedId], msa_3di[targetId], tmRes, map1, map2, qBtTM, tBtTM);
                float lddtTM  = getLDDTScore(seqDbrAA, seqDbr3Di, seqDbrCA, msaTM_aa,  par.pairThreshold);
                float lddt3di = getLDDTScore(seqDbrAA, seqDbr3Di, seqDbrCA, msa3di_aa, par.pairThreshold);
                msa_aa[mergedId]  = (lddtTM > lddt3di) ? msaTM_aa : msa3di_aa;
                msa_3di[mergedId] = (lddtTM > lddt3di) ? msaTM_3di : msa3di_3di;
                res               = (lddtTM > lddt3di) ? tmRes : res;
                qBt               = (lddtTM > lddt3di) ? qBtTM : qBt;
                tBt               = (lddtTM > lddt3di) ? tBtTM : tBt;
            } else {
                msa_aa[mergedId]  = msa3di_aa;
                msa_3di[mergedId] = msa3di_3di;
            }
            msa_aa[targetId] = "";
            msa_3di[targetId] = "";
            assert(msa_aa[mergedId].length() == msa_3di[mergedId].length());
            testSeqLens(msa_aa[mergedId], seqLens);
            updateBins(mergedId, targetId, res, map1, map2, qBt, tBt, neighbours); 

if (true) {
            // calculate LDDT of merged alignment
            float lddtScore = getLDDTScore(seqDbrAA, seqDbr3Di, seqDbrCA, msa_aa[mergedId], par.pairThreshold);
            std::cout << std::fixed << std::setprecision(3)
                << queryIsProfile << "\t" << targetIsProfile << '\t' << headers[mergedId] << "\t" << headers[targetId]
                << "\tLDDT: " << lddtScore << '\t' << res.score << '\n';
}

            // if (numMerges == 8) {
            //     exit(1);
            // }
            numMerges++;
            
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
                allSeqs_aa[mergedId]  = new Sequence(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_aa, 0, false, par.compBiasCorrection);
                allSeqs_3di[mergedId] = new Sequence(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, (const BaseMatrix *) &subMat_3di, 0, false, par.compBiasCorrection);
            }

            allSeqs_aa[mergedId]->mapSequence(mergedId, mergedId, profile_aa.c_str(), profile_aa.length() / Sequence::PROFILE_READIN_SIZE);
            allSeqs_3di[mergedId]->mapSequence(mergedId, mergedId, profile_3di.c_str(), profile_3di.length() / Sequence::PROFILE_READIN_SIZE);
            alreadyMerged[targetId] = true;
   
            // # Neighbours should == new sequence length
            // std::cout << neighbours[mergedId].size() << '\t' << allSeqs_3di[mergedId]->L << '\t' << allSeqs_aa[mergedId]->L << '\t' << mappings[mergedId].length() << '\n';
            // std::cout << profile_aa.length() << '\t' << Sequence::PROFILE_READIN_SIZE << '\t' << profile_aa.length() / Sequence::PROFILE_READIN_SIZE << '\n';
            // assert(neighbours[mergedId].size() == allSeqs_aa[mergedId]->L);
            // assert(neighbours[mergedId].size() == allSeqs_3di[mergedId]->L);
        }
        index += merges[i];
        // merged += merges[i];
    }

    // Find the final MSA (only non-empty string left in msa vectors)
    std::string finalMSA;
    for (size_t i = 0; i < sequenceCnt; ++i) {
        if (msa_aa[i] != "" && msa_3di[i] != "") {
            finalMSA = msa_aa[i];
            finalMSA_aa = msa_aa[i];
            finalMSA_3di = msa_3di[i];
            if (par.outputmode > 0) finalMSA = msa_3di[i];
            else finalMSA = msa_aa[i];
            continue;
        }
    }

    // Refine alignment -- MUSCLE5 style
    // 1. Partition into two sub-MSAs
    // 2. Remove all-gap columns
    // 3. Create sub-MSA profiles
    // 4. Save profiles -> Sequence objects
    // 5. Pairwise alignment
    // 6. Repeat x100
    if (par.refineIters > 0) {
        std::tie(finalMSA_aa, finalMSA_3di) = refineMany(
            tinySubMatAA,
            tinySubMat3Di,
            seqDbrAA,
            seqDbr3Di,
            seqDbrCA,
            finalMSA_aa,
            finalMSA_3di,
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
            par.evalProfile,
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
    }

    // Cleanup
    seqDbrAA.close();
    seqDbr3Di.close();
    delete[] alreadyMerged;
    delete [] tinySubMatAA;
    delete [] tinySubMat3Di;
    for (size_t i = 0; i < allSeqs_aa.size(); i++) {
        delete allSeqs_aa[i];
        delete allSeqs_3di[i];
    }
}

    // Write final MSA to file with correct headers
    DBWriter resultWriter(
        par.filenames[par.filenames.size()-1].c_str(),
        (par.filenames[par.filenames.size()-1] + ".index").c_str(),
        static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_OMIT_FILE
    );
    resultWriter.open();
    resultWriter.writeStart(0);
    std::string buffer;
    buffer.reserve(10 * 1024);
    KSeqWrapper* kseq = new KSeqBuffer(finalMSA_aa.c_str(), finalMSA_aa.length());
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
    FileUtil::remove((par.filenames[par.filenames.size()-1] + ".index").c_str());
   
    return EXIT_SUCCESS;
}

int structuremsa(int argc, const char **argv, const Command& command) {
    return structuremsa(argc, argv, command, false);
}

int structuremsacluster(int argc, const char **argv, const Command& command) {
    return structuremsa(argc, argv, command, true);
}
