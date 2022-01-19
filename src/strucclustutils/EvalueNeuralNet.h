//
// Created by Martin Steinegger on 19/01/2022.
//

#ifndef FOLDSEEK_EVALUENEURALNET_H
#define FOLDSEEK_EVALUENEURALNET_H
#include <cmath>
#include "kerasify/keras_model.h"
#include "BaseMatrix.h"

class EvalueNeuralNet {
private:
    BaseMatrix *subMat;
    double logDbResidueCount;
    KerasModel encoder;
    Tensor in;
    Tensor out;
public:

    EvalueNeuralNet(size_t dbResCount, BaseMatrix* subMat);

    std::pair<double, double> predictMuLambda(unsigned char * seq, unsigned int L);

    double computePvalue(double score, double lambda_, double mu) {
        double h = lambda_ * (score - mu);
        if(h > 10) {
            return -h;
        } else if (h < -2.5) {
            return -exp(-exp(-h));
        } else {
            return log((1.0 - exp(-exp(-h))));
        }
    }

    double computeEvalue(double score, double lambda_, double mu){
        return exp(computePvalue(score, lambda_, mu) + logDbResidueCount);
    }
};


#endif //FOLDSEEK_EVALUENEURALNET_H
