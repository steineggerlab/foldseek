#include "Matcher.h"
#include "MSANode.h"
#include <cassert>


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
    for (char c : str) if (match(c)) ++len;
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
 * Finds the index of an alignment start position in a pre-aligned sequence.
 * 
 * e.g. qStartPos = 2 in XXXX
 *
 *                    0 1   *
 *      Aligned = - - X X - X X
 *                0 1 2 3 4 5 6
 *
 *  5 in aligned string
 */
int findAlnStartIdx(std::string string, int start_pos) {
    int matches = 0;
    for (size_t i = 0; i < string.length(); ++i)
        if (!isGap(string[i]) && matches++ == start_pos)
            return i;
    return -1;
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
    return matches;
}

std::string noQueryMS(std::string targetSeq, int prevIdx, int curIdx) {
    std::cout << "\ntargetSeq: " << targetSeq;
    std::cout << "  prevIdx: " << prevIdx << std::endl;
    std::cout << "   curIdx: " << curIdx << std::endl;
    std::cout << "   length: " << targetSeq.length() << std::endl;

    assert(prevIdx >= 0);
    assert(curIdx >= prevIdx);
    assert(curIdx <= targetSeq.length());

    std::string extra;
    for (int i = prevIdx + 1; i <= curIdx; ++i) {
        std::cout << " char: " << targetSeq[i] << std::endl;
        if (!isGap(targetSeq[i]))
            extra.push_back(lwrchr(targetSeq[i]));
    }
    return extra;
}

std::string oneQueryMS(std::string targetSeq, int prevIdx, int curIdx) {
    // Add target sequence insertions, if any (lower cased)
    /* std::cout << "\ntargetSeq: " << targetSeq; */
    /* std::cout << "  prevIdx: " << prevIdx << std::endl; */
    /* std::cout << "   curIdx: " << curIdx << std::endl; */
    /* std::cout << "   length: " << targetSeq.length() << std::endl; */

    std::string extra;
    for (int i = prevIdx + 1; i < curIdx; ++i) {
        /* std::cout << " char: " << targetSeq[i] << std::endl; */
        if (!isGap(targetSeq[i]))
            extra.push_back(lwrchr(targetSeq[i]));
    }

    // Add match
    extra.push_back(uprchr(targetSeq[curIdx]));
    return extra;
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

// 
std::string alignNewSequence(std::string query, std::string target, size_t qId, size_t tId, size_t queryLen, Matcher::result_t result, std::vector<int> matches) {

    // TODO better way to initialise this?
    std::string newSeq = "";
    newSeq.reserve(32000);

    // Insert gaps up until query start position
    size_t j = 0;
    for (; j < (size_t)qId; ++j)
        newSeq.push_back('-');

    // Mark residues before target start position as insertion
    for (j = 0; j < (size_t)tId; ++j) {
        char c = target[j];
        c = (c == '-') ? '.' : lwrchr(c);
        newSeq.push_back(c);
    }

    // Check if alignment start in target is an insertion
    // TODO does this do anything? need a test case for it
    j = tId;
    while (insertion(target[j])) {
        int qIdx = matches[j];
        char residue = (target[j] == lwrchr(query[qIdx])) ? uprchr(target[j]) : lwrchr(target[j]);
        newSeq.push_back(residue);
        ++j;
    }
    tId = j;

    // Insert first match
    newSeq.push_back(target[tId]);

    std::cout << "  newSeq bef: " << newSeq << '\n';
    std::cout << "  matches: " ;
    for (int match : matches)
        std::cout << match << " ";
    std::cout << std::endl;
    
    int qPrev = qId;
    int tPrev = tId;

    for (int j = tId + 1; j <= result.dbEndPos; ++j) {

        // Get matching index in query
        int i = matches[j];
        assert(i >= 0);

        // Calculate offset between current and previous indices
        int dq = i - qPrev;
        int dt = j - tPrev;

        char matchState = result.backtrace[j - 1];
        std::cout
            << "  i: " << i << ", j: " << j << ", qPrev: " << qPrev << ", tPrev: " << tPrev << ", di: " << dq
            << ", dj: " << dt << ", qlen: " << query.length() << ", tlen: " << target.length()
            << ", state: " << matchState << "\n";
        assert(dq >= 0);
        assert(dt >= 0);

        if (dq == 1)
            newSeq.append(oneQueryMS(target, tPrev, j));
        else if (dq == 0)
            newSeq.append(noQueryMS(target, tPrev, j));
        else if (dq >= dt)
            gtQueryMS(newSeq, target, tPrev, j, dq, dt);
        else if (dt < dq)
            ltQueryMS(newSeq, target, tPrev, j, dq, dt);

        qPrev = i;
        tPrev = j;
        std::cout << "  newSeq aft: " << newSeq << '\n';
    }

    // Backfill gaps til equal match states
    size_t targetLen = countMatches(newSeq);
    assert(targetLen > 0);

    for (size_t j = targetLen; j < queryLen; ++j)
        newSeq.push_back('-');

    // add remaining seq. info as insertion state
    for (size_t j = tPrev + 1; j < target.length(); ++j) {
        if (target[j] == '\n')
            break;
        newSeq.push_back(target[j] == '-' ? '.' : lwrchr(target[j]));
    }

    // Save the new target to the current MSANode
    assert(target.length() > 0);

    // TODO maybe make new node here?
    newSeq.push_back('\n');
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
    /* summariseResult(result); */
    /* std::cout << " Query: " << query->sequence; */
    /* std::cout << "Target: " << target->sequence << std::endl; */
   
    // Identify index of first match state in query/target
    size_t qId = findAlnStartIdx(query->sequence, result.qStartPos);
    size_t tId = findAlnStartIdx(target->sequence, result.dbStartPos);

    // Map matching indices in query/target based on CIGAR
    // Need to be able to retrieve query index i given target index j when iterating through the backtrace
    std::vector<int> matches = mapBacktrace(query->sequence, target->sequence, qId, tId, result.backtrace);
    assert(matches.size() > 0);

    // Iterate all members of the new (target) MSANode, aligning them to the query (i.e. this MSANode)
    for (size_t i = 0; i < newNode->members.size(); ++i) {
        target = newNode->members[i];
        target->sequence = alignNewSequence(query->sequence, target->sequence, qId, tId, queryLen, result, matches);
        members.push_back(target);
    }

    return;
}

void MSANode::print() {
    for (MSASequence * node : members)
        node->print();
}
