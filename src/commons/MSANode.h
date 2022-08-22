#ifndef MSA_NODE_H
#define MSA_NODE_H

#include "Debug.h"
#include "Util.h"
#include "Matcher.h"

std::string noQueryMS(std::string targetSeq, int prevIdx, int curIdx);
std::string oneQueryMS(std::string targetSeq, int prevIdx, int curIdx);
void gtQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx, int dq, int dt);
void ltQueryMS(std::string &newSeq, std::string targetSeq, int prevIdx, int curIdx, int dq, int dt);

int findAlnStartIdx(std::string string, int start_pos);

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
    size_t queryLen,
    Matcher::result_t result,
    std::vector<int> matches
);

inline size_t countMatches(std::string str);

class MSASequence {
    public:
        size_t id;
        std::string sequence;

        MSASequence(
            size_t n_id,
            std::string n_sequence
        ) : id(n_id), sequence(n_sequence) {};
        MSASequence() = default;

        ~MSASequence() { };

        void print() {
            std::cout << '>' << id << '\n' << sequence;
        }
    private:
};

class MSANode {
    public:
        size_t id;
        std::vector<MSASequence *> members;

        MSANode(size_t c_id, std::vector<MSASequence *> c_members) : id(c_id), members(c_members) {
            alnLength = 0;
        }
        MSANode(size_t c_id) : id(c_id) {
            alnLength = 0;
        };
        MSANode() = default;

        MSASequence * findSeqById(size_t id);
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
