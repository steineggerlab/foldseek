#ifndef MSA_NODE_H
#define MSA_NODE_H

#include "Debug.h"
#include "Util.h"
#include "Matcher.h"

void noQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx);
void oneQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx);
void gtQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx, int dq, int dt);
void ltQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx, int dq, int dt);

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

int findNonGapIndex(std::string string, int pos);

std::vector<int> mapBacktrace(
    std::string query,
    std::string target,
    size_t q,
    size_t t,
    std::string backtrace
);

std::string alignNewSequence(
    std::string query,
    std::string target,
    size_t qId,
    size_t tId,
    size_t qEnd,
    size_t tEnd,
    size_t queryLen,
    Matcher::result_t result,
    std::vector<int> matches
);

inline size_t countMatches(std::string str);

class MSASequence {
    public:
        size_t id;
        int originalLen;
        std::string sequence;

        MSASequence(
            size_t n_id,
            std::string n_sequence
        ) : id(n_id), sequence(n_sequence) {
            originalLen = n_sequence.length();
        };
        MSASequence() = default;

        ~MSASequence() { };

        void print() {
            std::cout << '>' << id << '\n' << sequence;
        }
        std::string asString() {
            std::ostringstream oss;
            oss << '>' << id << '\n' << sequence << '\n';
            return oss.str();
        }
    private:
};

class MSANode {
    public:
        size_t id;
        int numMatches;
        std::vector<MSASequence *> members;

        MSANode(size_t c_id, std::vector<MSASequence *> c_members) : id(c_id), members(c_members) {
            alnLength = 0;
            numMatches = 0;
        };
        MSANode(size_t c_id) : id(c_id) {
            alnLength = 0;
            numMatches = 0;
        };
        MSANode() = default;

        MSASequence * findSeqById(size_t id);
        
        // Counts columns with % of gaps <= ratio 
        int countColumns(double ratio) {
            int count = 0;
            int gaps;
            size_t total = members.size();
            for (size_t i = 0; i < members[0]->sequence.length(); ++i) {
                gaps = 0;
                for (size_t j = 0; j < total; ++j)
                    gaps += (members[j]->sequence[i] == '-');
                if ((gaps / (double)total) <= ratio)
                    count++;
            }
            return count;
        }
        
        // Insert n gaps at residue id
        void insertGaps(size_t id, size_t n) {
            for (size_t i = 0; i < members.size(); ++i) {
                members[i]->sequence.insert(id, n, '-');
            }
        }

        void insertGapsAtEnd(size_t n) {
            for (size_t i = 0; i < members.size(); ++i) {
                members[i]->sequence.append(n, '-');
            }
        }
       
        void addNode(size_t merged_node_id, MSANode *new_cluster, Matcher::result_t &result);
        void print(); 
        void printMembers() {
            for (MSASequence *seq : members)
                std::cout << seq->id << '\n';
        };

    private:
        int alnLength;
};

#endif
