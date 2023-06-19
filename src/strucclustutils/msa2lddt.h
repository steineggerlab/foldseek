#include <vector>
#include <iostream>
#include "DBReader.h"
#include "KSeqWrapper.h"

void parseFasta(
    KSeqWrapper *kseq,
    DBReader<unsigned int> * seqDbrAA,
    DBReader<unsigned int> * seqDbr3Di,
    std::vector<std::string> &headers,
    std::vector<size_t>      &indices,
    std::vector<int>         &lengths,
    std::vector<std::string> &sequences,
    std::vector<std::string> &sequences3di,
    int &alnLength
);

std::tuple<std::vector<float>, std::vector<int>, float> calculate_lddt(
    std::vector<std::string> &sequences,
    std::vector<size_t> &indices,
    std::vector<int> &lengths,  // ungapped
    DBReader<unsigned int> * seqDbrCA,
    float pairThreshold
);

float getLDDTScore(
    DBReader<unsigned int> &seqDbrAA,
    DBReader<unsigned int> &seqDbr3Di,
    DBReader<unsigned int> &seqDbrCA,
    std::string msa,
    float pairThreshold
);