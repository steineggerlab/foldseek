#ifndef CALCPROBTP_H
#define CALCPROBTP_H

#include <cmath>

class CalcProbTP {
public:
  static float calculate(const float score){
  /*
   * Estimate probablity of a true positive (TP) given its score (structure bits)
   * See foldseek-analysis repo for the fitting procedure.
   */
    if (score <= 10){
      return 0;
    }
    if (score >= 100){
      return 1.0;
    }
    // Fitted score distributions for TPs and FPs
    float p_tp = (0.8279 * gamma_pdf(1.8123, 1/46.0042, score) + 0.1721 * gamma_pdf(1.0057, 1/563.5014, score)) * 0.1023;
    float p_fp = (0.34 * gamma_pdf(4.9259, 1/4.745, score) + 0.66 * gamma_pdf(9.4834, 1/1.3136, score)) * 0.8977;
    return 1 / (1 + (p_fp / p_tp));  // = p_tp / (p_tp + p_fp)
  }

private:
  static float gamma_pdf(const float alpha, const float beta, const float x){
    // Density of Gamma distribution
    // Parameters: alpha = shape, beta = 1 / scale
    return exp(alpha * log(beta) + (alpha-1) * log(x) + (-beta * x) - lgamma(alpha));
  }
};

#endif
