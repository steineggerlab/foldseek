#include "Matcher.h"
#include "MSANode.h"
#include <cassert>

#include <TMaligner.h>


inline bool match(char c) {
    return (c >= 'A' && c <= 'Z') || c == '-';
}

inline bool insertion(char c) {
    return (c >= 'a' && c <= 'z') || c == '.';
}

inline char uprchr(char chr) {
  return (chr >= 'a' && chr <= 'z') ? chr + 'A' - 'a' : chr;
}

inline char lwrchr(char chr) {
    return (chr >= 'A' && chr <= 'Z') ? chr - 'A' + 'a' : chr;
}

inline char isGap(char chr) {
    return (chr == '-' || chr == '.');
}

// Find index of next non-gap character in alignment string
inline size_t findNextResidueIdx(std::string str, size_t start) {
    ++start;
    while (start < str.length() && isGap(str[start]))
        ++start;
    return start;
}

// Count match states in an aligned sequence
inline size_t countMatches(std::string str) {
    size_t len = 0;
    for (char c : str)
        if (match(c))
            ++len;
    return len;
}

void summariseResult(Matcher::result_t result) {
    std::cout << "Result:\n";
    std::cout << "        dbKey: " << result.dbKey << std::endl;
    std::cout << "        Score: " << result.seqId * 100 << " " << result.eval << std::endl;
    std::cout << "  qStart/qEnd: " << result.qStartPos << " " << result.qEndPos << std::endl;
    std::cout << "  tStart/tEnd: " << result.dbStartPos << " " << result.dbEndPos << std::endl;
    std::cout << "        CIGAR: " << result.backtrace << " (" << result.backtrace.length() << ")\n" << std::endl;
}

/**
 * Finds residue index from unaligned sequence in pre-aligned sequence.
 * 
 * e.g. qStartPos = 2 in XXXX
 *
 *                    0 1   *
 *      Aligned = - - X X - X X
 *                0 1 2 3 4 5 6
 *
 *  5 in aligned string
 */
int findNonGapIndex(std::string string, int pos) {
    int index = -1;
    int matches = 0;
    for (size_t i = 0; i < string.length(); ++i) {
        if (!isGap(string[i])) {
            if (matches == pos) {  
                index = matches;
                break;
            }
            ++matches;
        } 
    }
    return index;
}

void print_vector(std::vector<int> vector) {
    for (int member : vector)
        std::cout << member << " ";
    std::cout << std::endl;
}

// Map matching indices in query/target based on CIGAR
// Need to be able to retrieve query index i given target index j when iterating through the backtrace
std::vector<int> mapBacktrace(std::string query, std::string target, size_t q, size_t t, std::string backtrace) {
    std::vector<int> matches(target.length());
    assert(matches.size() > 0);
    size_t b = 0;
    size_t qq = q;      // NOTE: keeping track of the old q here works..
    while (b < backtrace.length()) {
        if (t >= target.length() || q >= query.length())
            break;
        switch (backtrace[b]) {
            case 'M': {
                matches[t] = qq = q;
                size_t newT = findNextResidueIdx(target, t);
                for (size_t x = t; x < newT; ++x)
                    matches[x] = qq;
                t = newT;
                q = findNextResidueIdx(query, q);
                break;
            }
            case 'I': {
                qq = q;
                q = findNextResidueIdx(query, q);
                break;
            }
            case 'D': {
                size_t newT = findNextResidueIdx(target, t);
                for (size_t x = t; x < newT; ++x)
                    matches[x] = qq;
                t = newT;
                break;
            }
        }
        ++b;
        assert(t <= matches.size());
    }

    // TODO maybe not necessary
    // Fill in remainder of indices
    
    for (; t < target.length(); ++t)
        matches[t] = q;

    return matches;
}

void noQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx) {
    assert(prevIdx >= 0);
    assert(curIdx >= prevIdx);
    assert(curIdx <= targetSeq.length());

    for (int i = prevIdx + 1; i <= curIdx; ++i) {
        if (!isGap(targetSeq[i]))
            newSeq.push_back(lwrchr(targetSeq[i]));
    }
}

void oneQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx) {
    // Add target sequence insertions, if any (lower cased)
    /* std::cout << "\ntargetSeq: " << targetSeq; */
    /* std::cout << "  prevIdx: " << prevIdx << std::endl; */
    /* std::cout << "   curIdx: " << curIdx << std::endl; */
    /* std::cout << "   length: " << targetSeq.length() << std::endl; */
    for (int i = prevIdx + 1; i < curIdx; ++i) {
        /* std::cout << " char: " << targetSeq[i] << std::endl; */
        if (!isGap(targetSeq[i]))
            newSeq.push_back(lwrchr(targetSeq[i]));
    }
    // Add match
    newSeq.push_back(targetSeq[curIdx]);
}

void gtQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx, int dq, int dt) {
    // Calculate bound 
    int bound = prevIdx + (int)(dt / 2);

    // Add left-bounded residues
    for (int i = prevIdx + 1; i <= bound; ++i)
        newSeq.push_back(uprchr(targetSeq[i]));

    // Add central gaps
    for (int gap = 1; gap <= dq - dt; ++gap)
        newSeq.push_back('-');

    // Add right-bounded residues
    for (int i = bound + 1; i <= curIdx; ++i){
        newSeq.push_back(uprchr(targetSeq[i]));
    }
}

void ltQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx, int dq, int dt) {
    int bound = prevIdx + (int)(dt / 2);

    // Add left-bounded template residues
    for (int i = prevIdx + 1; i < bound; ++i)
        newSeq.push_back(uprchr(targetSeq[i]));

    // Add central inserts
    for (int gap = 1; gap <= dt - dq; gap++, bound++) {
        char c = targetSeq[bound];
        if (!isGap(c))
            newSeq.push_back(lwrchr(c));
    }

    // Add right-bounded residues
    for (int i = bound; i <= curIdx; ++i)
        newSeq.push_back(uprchr(targetSeq[i]));
}

std::string alignWithTMAlign(std::string query, std::string target, Matcher::result_t result) {
    // TMaligner::TMaligner aligner;

    // Pre: take only matching residues in query/target based on Foldseek backtrace
    // e.g. ABCDEFGHI     ABCDGHI
    //      MMMMIIMMM --> MMMMMMM
    //      ABCD--GHI     ABCDGHI
    int q = findNonGapIndex(query, result.qStartPos);
    int t = findNonGapIndex(target, result.dbStartPos);

    size_t numMatches = (size_t)std::count(result.backtrace.begin(), result.backtrace.end(), 'M');
    char newQuery[numMatches];
    char newTarget[numMatches];

    size_t b = 0;
    size_t m = 0;
    for (; m < numMatches && b < result.backtrace.length(); ++b) {
        switch (result.backtrace[b]) {
            case 'M': {
                newQuery[m] = query[q];
                newTarget[m] = target[t];
                ++m;
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
    
    std::cout << "     Query: " << query;
    std::cout << "    Target: " << target;
    std::cout << " New query: " << newQuery << std::endl;
    std::cout << "New target: " << newTarget << std::endl;
    
    // Initialise query on the aligner (xyz, sequence, length)
    // aligner.initQuery();
    
    // Align using TMAlign
    // Matcher::result_t result = aligner.align();
    
    // Adjust target based on new backtrace

    return "abc";
}

// 
std::string alignNewSequence(std::string query, std::string target, size_t qId, size_t tId, size_t qEnd, size_t tEnd, size_t queryLen, Matcher::result_t result, std::vector<int> matches) {
    // TODO better way to initialise this?
    std::string newSeq = "";
    newSeq.reserve(32000);

    // Insert gaps up until query start position
    size_t j = 0;
    for (; j < qId; ++j)
        newSeq.push_back('-');
    
    // TODO might have to calculate this for each target sequence
    // tId = findAlnStartIdx(target, result.dbStartPos);

    for (j = 0; j < tId; ++j) {
        char c = target[j];
        c = (c == '-') ? '-' : lwrchr(c);
        newSeq.push_back(c);
    }

    // TODO does this do anything? need a test case for it
    //      Since I already calculate start pos accounting for dels, maybe not
    // Check if alignment start in target is an insertion
    j = tId;
    while (insertion(target[j])) {
        int qIdx = matches[j];
        char residue = (target[j] == lwrchr(query[qIdx])) ? uprchr(target[j]) : lwrchr(target[j]);
        newSeq.push_back(residue);
        ++j;
    }
    tId = j;

    // Insert first match
    newSeq.push_back((islower(target[tId])) ? uprchr(target[tId]) : target[tId]);

    // std::cout << "  newSeq bef: " << newSeq << '\n';
    std::cout << "  matches: " ;
    for (int match : matches)
        std::cout << match << " ";
    // std::cout << std::endl;
    // std::cout << "current sequence: " << target << std::endl;
    
    int qPrev = qId;
    int tPrev = tId;

    // TODO dbEndPos relates only to the target sequence driving this merge
    // and may be longer than the current sequence being merged
    // - all sequences in a node should be same size, so should fix by ensuring
    //   alignments are same length at end of alignment
    int b = 1;
    for (size_t j = tPrev + 1; j <= tEnd; ++j) {

        // Get matching index in query
        int i = matches[j];
        assert(i >= 0);

        // Calculate offset between current and previous indices
        int dq = i - qPrev;
        int dt = j - tPrev;

        char matchState = result.backtrace[b++];
        std::cout
            << "  i: " << i << ", j: " << j << ", qPrev: " << qPrev << ", tPrev: " << tPrev << ", di: " << dq
            << ", dj: " << dt << ", qlen: " << query.length() << ", tlen: " << target.length()
            << ", state: " << matchState << "\n";
        assert(dq >= 0);
        assert(dt >= 0);

        if (dq == 1)
            oneQueryMS(newSeq, target, tPrev, j);
        else if (dq == 0)
            noQueryMS(newSeq, target, tPrev, j);
        else if (dq >= dt)
            gtQueryMS(newSeq, target, tPrev, j, dq, dt);
        else if (dt < dq)
            ltQueryMS(newSeq, target, tPrev, j, dq, dt);

        qPrev = i;
        tPrev = j;
        // std::cout << "  newSeq aft: " << newSeq << '\n';
    }
   
    std::string::difference_type numMatches = std::count(result.backtrace.begin(), result.backtrace.end(), 'M');
    
    // Backfill gaps til equal match states
    size_t targetLen = countMatches(newSeq);
    assert(targetLen > 0);
    std::cout<< "queryLen: " << queryLen << ", targetLen: " << targetLen << ", numMatches: " << numMatches << std::endl;
    
    int numGaps = queryLen - result.qEndPos - 1;        // How many gaps to insert to match query
    std::cout << ", numGaps: " << numGaps << ", qEndPos: " << result.qEndPos + 1 << std::endl;

    for (size_t j = 0; j < numGaps; ++j)
        newSeq.push_back('-');

    // for (size_t j = targetLen; j < (size_t)numMatches; ++j)
    //     newSeq.push_back('-');
    
    // add remaining seq. info as insertion state
    for (size_t j = tPrev + 1; j < target.length(); ++j) {
        newSeq.push_back(target[j]);
    }

    // for (size_t j = tPrev + 1; j < target.length(); ++j) {
    //     if (target[j] == '\n')
    //         break;
    //     newSeq.push_back(target[j] == '-' ? '-' : target[j]);
    // }
    
    // Add the remaining gaps '-' to the end of the template sequence
    // for (i = hit.i2 + 1; i <= L; ++i) {
    //     cur_seq[h++] = '-';
    //     // too few columns? Reserve double space
    //     if (h >= maxcol - 1000) {
    //         char *new_seq = new char[2 * maxcol];
    //         strncpy(new_seq, cur_seq, h);
    //         delete[] (cur_seq);
    //         cur_seq = new_seq;
    //         maxcol *= 2;
    //     }
    // }
    // // add remaining seq. info as insertion state
    // while(Tali.seq[k][ll] != '\0'){
    //     if(Tali.seq[k][ll] == '-'){
    //         cur_seq[h++] = '.';
    //     }else{
    //         cur_seq[h++] = lwrchr(Tali.seq[k][ll]);
    //     }
    //     ll++;
    // }
    // cur_seq[h++] = '\0';

    // Save the new target to the current MSANode
    assert(target.length() > 0);

    // TODO index is likely pushing too far forward, copying \n from original target sequence AND adding this new one
    // newSeq.push_back('\n');
    
    // std::cout << "final newSeq: " << newSeq;

    return newSeq;
}

MSASequence * MSANode::findSeqById(size_t seqId) {
    for (size_t i = 0; i < members.size(); ++i)
        if (members[i]->id == seqId)
            return members[i];

    Debug(Debug::ERROR) << "Failed to find node with id " << seqId << '\n';
    EXIT(EXIT_FAILURE);
}

/**
 * Merges a new cluster into the current cluster using a given result. 
 *
 */
void MSANode::addNode(size_t queryID, MSANode *newNode, Matcher::result_t &result) {
    /* std::cout << "-----------------\nStart add cluster\n-----------------\n"; */

    // Grab MSASequence objects corresponding to query/target in result
    MSASequence *query = findSeqById(queryID);
    MSASequence *target = newNode->findSeqById(result.dbKey);

    // Count number of match states in query sequence
    // Used to calculate gaps at end of target
    size_t queryLen = countMatches(query->sequence);
    assert(queryLen > 0);
    if (alnLength == 0)
        alnLength = queryLen;

    // Result summary
    // summariseResult(result);
    std::cout << " Query " << query->id << ": " << query->sequence << " " << query->sequence.length() << std::endl;
    std::cout << "Target " << target->id << ": " << target->sequence << " " << target->sequence.length() << std::endl;
   
    // Identify index of first match state in query/target
    size_t qId = findNonGapIndex(query->sequence, result.qStartPos);
    size_t tId = findNonGapIndex(target->sequence, result.dbStartPos);
    size_t qEnd = findNonGapIndex(query->sequence, result.qEndPos);
    size_t tEnd = findNonGapIndex(target->sequence, result.dbEndPos);
    std::cout << result.qStartPos << " " << result.qEndPos << std::endl;
    std::cout << qId << " " << qEnd << std::endl;
    std::cout << result.dbStartPos << " " << result.dbEndPos << std::endl;
    std::cout << tId << " " << tEnd << std::endl;
    std::cout << "CIGAR: " << result.backtrace << std::endl;

    // Map matching indices in query/target based on CIGAR
    // Need to be able to retrieve query index i given target index j when iterating through the backtrace
    std::vector<int> matches = mapBacktrace(query->sequence, target->sequence, qId, tId, result.backtrace);
    assert(matches.size() > 0);

    // Iterate all members of the new (target) MSANode, aligning them to the query (i.e. this MSANode)
    for (size_t i = 0; i < newNode->members.size(); ++i) {
        target = newNode->members[i];
        target->sequence = alignNewSequence(query->sequence, target->sequence, qId, tId, qEnd, tEnd, queryLen, result, matches);
        // alignWithTMAlign(query->sequence, target->sequence, result);
        members.push_back(target);
    }

    return;
}

void MSANode::print() {
    for (MSASequence * seq : members) {
        seq->print();
        std::cout << std::endl;
    }
}
