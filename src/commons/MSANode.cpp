#include "Matcher.h"
#include "MSANode.h"
#include "filesystem"
#include <cassert>

#include <TMaligner.h>

// Map residue in query to target
struct ResidueMapping {
    char state;
    int qPos;
    int tPos;
    bool insertion;
    
    ResidueMapping(char state, int qPos, int tPos, bool insertion) : state(state), qPos(qPos), tPos(tPos), insertion(insertion) {};
    ResidueMapping() = default;
    ~ResidueMapping() {};
    
    void print() {
        std::cout << "State: " << state << ", qPos: " << qPos << ", tPos: " << tPos << ", Insertion: " << insertion << std::endl;
    }
};

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
                index = i;
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

/**
 * @brief Aligns two MSANodes using ResidueMappings from Foldseek backtraces
 * 
 * @param query - query MSANode*
 * @param target - target MSANode*
 * @param mapping - mappings between query and target residue indices
 */
void alignNodes(MSANode *query, MSANode *target, std::vector<ResidueMapping> mapping) {
    
    ResidueMapping curRM;
    ResidueMapping oldRM;
    
    // Keep track of new gaps in query/target to adjust indices
    int qGaps = mapping[0].tPos;
    int tGaps = mapping[0].qPos;
    
    // Also characters after final match position before sequences are modified
    int tEndChars = target->members[0]->sequence.length() - 1 - mapping[mapping.size() - 1].tPos;
    int qEndChars = query->members[0]->sequence.length() - 1 - mapping[mapping.size() - 1].qPos;

    // Insert start gaps
    // Query: qPos gaps from 0
    // Target: tPos gaps from qPos
    //               v 
    // e.g. Q: ---xxxXXXXX
    //      T: yyy---YYYYY
    query->insertGaps(0, qGaps);
    target->insertGaps(qGaps, tGaps); 
    
    for (size_t i = 1; i < mapping.size(); ++i) {
        curRM = mapping[i];
        oldRM = mapping[i - 1];
        
        int dq = curRM.qPos - oldRM.qPos;
        int dt = curRM.tPos - oldRM.tPos;
        
        if (dq == 1) {
            // One match in query
            // Mark dt - 1 characters before match as insertion
            query->insertGaps(oldRM.qPos + qGaps + 1, dt - 1);
            qGaps += dt - 1;

        } else if (dq == 0) {
            // No matches in query, all target residues are insertions
            // Insert dt gaps in query
            query->insertGaps(oldRM.qPos + qGaps + 1, dt);
            qGaps += dt;

        } else if (dq >= dt) {
            // More query matches than target matches
            // Insert dq - dt gaps in target
            // TODO insert at middle of oldRM.tPos and curRM.tPos
            target->insertGaps(oldRM.tPos + tGaps + 1, dq - dt);
            tGaps += dq - dt;
            
        } else if (dt >= dq) {
            // More target matches than query matches
            // Insert dt - dq gaps in query
            // TODO mark from middle of oldRM.tPos and curRM.tPos
            query->insertGaps(oldRM.qPos + qGaps + 1, dt - dq);
            qGaps = dt - dq;
        }
    }
    
    // End gaps
    // Query: from dbEndPos to end of target
    // Target: from qEndPos to end of query
    if (tEndChars > 0)
        query->insertGapsAtEnd(tEndChars);
    if (qEndChars > 0)
        target->insertGaps(curRM.tPos + 1 + tGaps, qEndChars);
   
    // Add target members to query
    query->members.insert(query->members.end(), target->members.begin(), target->members.end());
}

// Map matching indices in query/target based on CIGAR
// Need to be able to retrieve query index i given target index j when iterating through the backtrace
std::vector<ResidueMapping> mapBacktraceRM(std::string query, std::string target, size_t q, size_t t, std::string backtrace) {
    std::vector<ResidueMapping> matchesRM;
    matchesRM.reserve(backtrace.length());
    size_t b = 0;
    while (b < backtrace.length()) {
        if (t >= target.length() || q >= query.length())
            break;
        switch (backtrace[b]) {
            case 'M': {
                matchesRM.emplace_back('M', q, t, islower(query[q]));
                t = findNextResidueIdx(target, t);
                q = findNextResidueIdx(query, q);
                break;
            }
            case 'I': {
                q = findNextResidueIdx(query, q);
                break;
            }
            case 'D': {
                t = findNextResidueIdx(target, t);
                break;
            }
        }
        ++b;
    }
    return matchesRM;
}

// Map matching indices in query/target based on CIGAR
// Need to be able to retrieve query index i given target index j when iterating through the backtrace
std::vector<int> mapBacktrace(std::string query, std::string target, size_t q, size_t t, std::string backtrace) {
    std::cout << "Mapping backtrace\n"; 
    std::cout << "   Query: " << query << std::endl;
    std::cout << "  Target: " << target << std::endl;
    std::cout << "       q: " << q << std::endl;
    std::cout << "       t: " << t << std::endl;

    std::vector<int> matches(target.length());

    assert(matches.size() > 0);
    size_t b = 0;
    size_t qq = q;
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
    // std::cout << "          Adding: " << uprchr(targetSeq[curIdx]) <<std::endl;
    newSeq.push_back(uprchr(targetSeq[curIdx]));
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

std::string alignNewSequence(std::string query, std::string target, size_t qId, size_t tId, size_t qEnd, size_t tEnd, size_t queryLen, Matcher::result_t result, std::vector<int> matches) {

    // TODO better way to initialise this?
    std::string newSeq = "";
    newSeq.reserve(query.length());

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
    // std::cout << "  matches: " ;
    // for (int match : matches)
    //     std::cout << match << " ";
    // std::cout << std::endl;
    // std::cout << "current sequence: " << target << std::endl;
    
    int qPrev = qId;
    int tPrev = tId;
    
    // TODO insertions getting marked single index too early

    // TODO dbEndPos relates only to the target sequence driving this merge
    // and may be longer than the current sequence being merged
    // - all sequences in a node should be same size, so should fix by ensuring
    //   alignments are same length at end of alignment
    int b = 1;
    for (size_t j = tPrev + 1; j <= tEnd; ++j) {
        if (j >= target.length()) {
            newSeq += '-';
            int i = matches[j];
            int dq = i - qPrev;
            int dt = j - tPrev;
            qPrev = i;
            tPrev = j;
            continue;
        }

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
            << ", state: " << matchState;
        assert(dq >= 0);
        assert(dt >= 0);

        if (dq == 1) {
            std::cout << ", oneQueryMS";
            oneQueryMS(newSeq, target, tPrev, j);
        } else if (dq == 0) {
            std::cout << ", noQueryMS";
            noQueryMS(newSeq, target, tPrev, j);
        } else if (dq >= dt) {
            std::cout << ", gtQueryMS";
            gtQueryMS(newSeq, target, tPrev, j, dq, dt);
        } else if (dt < dq) {
            std::cout << ", ltQueryMS";
            ltQueryMS(newSeq, target, tPrev, j, dq, dt);
        }
        std::cout << ", new char: " << newSeq[newSeq.length() - 1] << std::endl;

        qPrev = i;
        tPrev = j;
        // std::cout << "  newSeq aft: " << newSeq << '\n';
    }
    
  
    std::string::difference_type numMatches = std::count(result.backtrace.begin(), result.backtrace.end(), 'M');
    
    // Backfill gaps til equal match states
    size_t targetLen = countMatches(newSeq);
    assert(targetLen > 0);
    std::cout << "\nqueryLen: " << queryLen << ", targetLen: " << targetLen << ", numMatches: " << numMatches;
    
    int numGaps;
    if (queryLen > targetLen)
        numGaps = queryLen - targetLen;        // How many gaps to insert to match query
    else
        numGaps = 0;

    std::cout << ", numGaps: " << numGaps << ", qEndPos: " << result.qEndPos + 1 << std::endl;
    std::cout << std::endl;
    std::cout << newSeq << std::endl;
    
    // Gaps to equal query match states
    for (size_t j = 0; j < numGaps; ++j)
        newSeq.push_back('-');

    // Add remaining sequence
    for (size_t j = tPrev + 1; j < target.length(); ++j) {
        newSeq.push_back(lwrchr(target[j]));
    }

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
    numMatches = queryLen;
    assert(queryLen > 0);
    if (alnLength == 0)
        alnLength = queryLen;
    
    // Identify index of first match state in query/target
    size_t qId = findNonGapIndex(query->sequence, result.qStartPos);
    size_t tId = findNonGapIndex(target->sequence, result.dbStartPos);
    size_t qEnd = findNonGapIndex(query->sequence, result.qEndPos);
    size_t tEnd = findNonGapIndex(target->sequence, result.dbEndPos);

    // Map matching indices in query/target based on CIGAR
    // std::vector<int> matches = mapBacktrace(query->sequence, target->sequence, qId, tId, result.backtrace);
    std::vector<ResidueMapping> matchesRM = mapBacktraceRM(
        query->sequence,
        target->sequence,
        qId,
        tId,
        result.backtrace
    );
    assert(matchesRM.size() > 0);
    
    // Iterate all members of the new (target) MSANode, aligning them to the query (i.e. this MSANode)
    // for (size_t i = 0; i < newNode->members.size(); ++i) {
    //     target = newNode->members[i];
    //     target->sequence = alignNewSequence(query->sequence, target->sequence, qId, tId, qEnd, tEnd, queryLen, result, matches);
    //     target->print();
    //     members.push_back(target);
    // }
    
    alignNodes(this, newNode, matchesRM);
    
    // Run reformat.pl to make sure a3m is valid after every merge
    // {
    //     char fileName[] = "/tmp/mytemp.XXXXXX";
    //     int fd = mkstemp(fileName);
        
    //     std::stringstream ss;
    //     for (MSASequence *seq : members)
    //         ss << '>' << seq->id << '\n' << seq->sequence << '\n';
    //     std::string alnFASTA = ss.str();
    //     std::cout << "\n#######\n";
    //     std::cout << alnFASTA << std::endl;
    //     std::cout << "\n#######\n";
    //     write(fd, alnFASTA.c_str(), alnFASTA.size());
        
    //     char fileNameOut[] = "/tmp/mytemp.XXXXXX";
    //     int fdOut = mkstemp(fileNameOut);  // Necessary?
        
    //     std::string command = std::string("perl ~/Downloads/reformat.pl a3m fas ") + fileName + " " + fileNameOut;
    //     int code = system(command.c_str());
    //     unlink(fileName);
    //     unlink(fileNameOut);
    //     if (code == -1 || WEXITSTATUS(code) != 0) {
    //         std::cout << "Error when reformatting alignment" << std::endl;
    //         exit(code);
    //     }
    //     std::cout << "return code after reformat.pl: " << code << std::endl;
    // }
    
    return;
}

void MSANode::print() {
    for (MSASequence * seq : members) {
        seq->print();
        std::cout << std::endl;
    }
}
