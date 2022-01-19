//
// Created by Martin Steinegger on 19/01/2022.
//

#include "EvalueNeuralNet.h"
#include "evalue_nn.kerasify.h"


EvalueNeuralNet::EvalueNeuralNet(size_t dbResCount, BaseMatrix* subMat) : subMat(subMat) {
        logDbResidueCount = log(static_cast<double>(dbResCount));
        encoder.LoadModel(
        std::string((const char *)evalue_nn_kerasify,
        evalue_nn_kerasify_len));
        in = Tensor(subMat->alphabetSize + 1);
        out = Tensor(2);
}

std::pair<double, double> EvalueNeuralNet::predictMuLambda(unsigned char * seq, unsigned int L){
    for(int i = 0; i < subMat->alphabetSize; i++){
        in.data_[i] = 0;
    }
    for (unsigned int i = 0; i < L; i++) {
        in.data_[seq[i]]++; ;
    }
    in.data_[subMat->alphabetSize] = L;
    encoder.Apply(&in, &out);
    // used to normalize the output
    double mu1 = 0.25319026504505765;
    double sigma1 = 0.0870925635266098;
    double mu2 = 16.48226262676655;
    double sigmal2 = 3.472243709346458;
    return std::make_pair(out.data_[0]*sigma1+mu1,
                          out.data_[1]*sigmal2+mu2);
}