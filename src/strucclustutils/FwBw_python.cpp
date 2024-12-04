#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include <cstring>
#include <sstream>
#include <chrono>

extern "C" {
void align(double* scoreForward, double* P, double* S_log,
            size_t queryLen, size_t targetLen,
            float go, float ge, float T, int length, int blocks);
double** allocateMemory(size_t queryLen, size_t targetLen);
void rescaleBlocks(double **matrix, double **scale, size_t rows, size_t length, size_t blocks, size_t targetLen);
void forwardBackwardSaveBlockMaxLocal(double* S,double* S_log, double** z_init,
                                                   float T, float go, float ge,
                                                   size_t rows, size_t start, size_t end, size_t memcpy_cols, size_t targetlen,
                                                   double** zm, double* zmax, double** zmBlock, double* zeBlock, double* zfBlock);

}


void align(double* scoreForward, double* P, double* S_log,
            size_t queryLen, size_t targetLen,
            float go, float ge, float T, int length, int blocks) {


    
    double* scoreBackward = new double[queryLen * targetLen];
    double** zmForward = allocateMemory(queryLen, targetLen);
    double** zmBackward = allocateMemory(queryLen, targetLen);
    double** zmBlock = allocateMemory(queryLen + 1, length + 1);
    double** zmaxForward = allocateMemory(blocks, queryLen);
    double** zmaxBackward = allocateMemory(blocks, queryLen);
    double* zeBlock = new double[length + 1];
    double* zfBlock = new double[length + 1];


    for(size_t i = 0; i < queryLen; ++i){
        for(size_t j = 0; j < targetLen; ++j){
            scoreBackward[i * targetLen + j] = scoreForward[(queryLen - 1 - i)*targetLen + (targetLen - 1 - j)];
        }
    }


   /*for (size_t j = 0; j < targetLen; ++j) {
            zmForward[0][j] = -DBL_MAX;
            zmBackward[0][j] = -DBL_MAX;

    }
    for(size_t i = 0; i < queryLen; ++i){
        zmForward[i][0] = -DBL_MAX;
        zmBackward[i][0] = -DBL_MAX;
    }*/


    double* zInit[3];
    zInit[0] = new double[queryLen];
    zInit[1] = new double[queryLen];
    zInit[2] = new double[queryLen];


    for (unsigned int i=0; i < queryLen; ++i){
        zInit[0][i] = -DBL_MAX;
        zInit[1][i] = -DBL_MAX;
        zInit[2][i] = -DBL_MAX;
    }
 
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;
        forwardBackwardSaveBlockMaxLocal(scoreForward, S_log, zInit, T, go, ge, queryLen, start, end, memcpy_cols, targetLen,
                                         zmForward, zmaxForward[b], zmBlock ,zeBlock, zfBlock);
        
    }
    ///////////////////////////////////Backward////////////////////////////////////////
    for (unsigned int i=0; i < queryLen; ++i){
        zInit[0][i] = -DBL_MAX;
        zInit[1][i] = -DBL_MAX;
        zInit[2][i] = -DBL_MAX;
    }

    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, static_cast<size_t>(targetLen)) - start;
        
        forwardBackwardSaveBlockMaxLocal(scoreBackward, S_log, zInit, T, go, ge, queryLen, start, end, memcpy_cols, targetLen,
                                         zmBackward, zmaxBackward[b], zmBlock, zeBlock, zfBlock);
    }

    ///////////////////////////////////Rescale////////////////////////////////////////
    // Rescale the values by the maximum in the log space for each block
    // This turns the matrix into log space

    
    rescaleBlocks(zmForward, zmaxForward, queryLen, length, blocks, targetLen);
    rescaleBlocks(zmBackward, zmaxBackward, queryLen, length, blocks, targetLen);
    /*std::ofstream outfile;
    outfile.open("/home/lasse/Desktop/Projects/FB_martin/zmForward.txt");
    if (outfile.is_open()) {
        for (size_t i = 0; i < queryLen; ++i) {
            for (size_t j = 0; j < targetLen; ++j) {
                outfile << zmForward[i][j] << " ";
            }
            outfile << "\n";
        }
        outfile.close();
    } else {
        std::cout << "Unable to open file";
    }*/


    double max_zm = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            max_zm = std::max(max_zm, zmForward[i][j]);
        }
    }

    double sum_exp= 0.0;
    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            sum_exp += exp(zmForward[i][j] - max_zm);
        }
    }
    double logsumexp_zm = max_zm + log(sum_exp);




    // compute posterior probabilities

    //float maxP = -std::numeric_limits<float>::max();
    
    //outfile.open("/home/lasse/Desktop/Projects/FB_martin/hello.txt");


    for (size_t i = 0; i < queryLen; ++i) {
        for (size_t j = 0; j < targetLen; ++j) {
            // Debug(Debug::INFO) << zmForward[i][j] << '\t' << zmBackward[queryLen - 1 - i][targetLen - 1 - j] << '\n';
            P[i*targetLen + j] = exp(
                zmForward[i][j]
                + zmBackward[queryLen - 1 - i][targetLen - 1 - j]
                - S_log[i*targetLen + j] // FIXME scoreForward is already exp(S/T)
                - logsumexp_zm
            );
            
            //outfile << i << "," << j << "," <<  zmForward[i][j] << "," << zmBackward[queryLen - 1 - i][targetLen - 1 - j] << std::endl;
            
            //maxP = std::max(maxP, P[i*targetLen + j]);
            // Debug(Debug::INFO) << P[i][j] << '\t';
        }
        // Debug(Debug::INFO) << '\n';
    }
    //outfile.close();


    free(scoreBackward);
    free(zmForward);
    free(zmBackward);
    free(zmBlock);
    free(zmaxForward);
    free(zmaxBackward);
    free(zeBlock);
    free(zfBlock);
    return;
}


void forwardBackwardSaveBlockMaxLocal(double* S, double* S_log, double** z_init,
                                                   float T, float go, float ge,
                                                   size_t rows, size_t start, size_t end, size_t memcpy_cols, size_t targetlen,
                                                   double** zm, double* zmax, double** zmBlock, double* zeBlock, double* zfBlock) {
    double exp_go = exp(go / T);
    double exp_ge = exp(ge / T);
    
    
    

    //std::cout << "test" << std::endl;
    memset(zeBlock, 0, (end - start + 1) * sizeof(double)); 
    memset(zfBlock, 0, (end - start + 1) * sizeof(double)); 
 
    std::vector<double> ze_first(rows+1, 0);
    std::vector<double> zf_first(rows+1, 0);
    
    //Init blocks
    memset(zmBlock[0], 0, (end - start + 1) * sizeof(double));

    for (size_t i = 0; i < rows; ++i) {
        zmBlock[i+1][0] = z_init[0][i];
        ze_first[i+1] = z_init[1][i];
        zf_first[i+1] = z_init[2][i];
    }


    size_t cols = memcpy_cols;


    double current_max = 0;
    for (size_t i = 1; i <= rows; ++i) {
        if (i != 1) {
            zmBlock[i - 1][0] = exp(zmBlock[i - 1][0]);
            ze_first[i - 1] = exp(ze_first[i - 1]);
            zf_first[i - 1] = exp(zf_first[i - 1]);
            // Debug(Debug::INFO) << zmBlock[i - 1][0] << '\t';
        }
        const float expMax = exp(-current_max);
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            if (j == 1){
                double tmp = (zmBlock[i - 1][j - 1] + ze_first[i-1] + zf_first[i - 1] + expMax);
                zmBlock[i][j] = tmp * S[(i - 1) * targetlen + start + j - 1];
            }
            else{
                double tmp = (zmBlock[i - 1][j - 1] + zeBlock[j - 1] + zfBlock[j - 1] + expMax);
                zmBlock[i][j] = tmp * S[(i - 1) * targetlen + start + j - 1];
            }
        }
        

        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            if (j == 1) {
                zeBlock[j] = exp(zmBlock[i][j - 1]) * exp_go + exp(ze_first[i]) * exp_ge;
            } else {
                zeBlock[j] = zmBlock[i][j - 1] * exp_go + zeBlock[j - 1] * exp_ge;
            }

        }
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
                zfBlock[j] = zmBlock[i - 1][j] * exp_go + zfBlock[j] * exp_ge;
        }

        double z_temp = *std::max_element(zmBlock[i] + 1, zmBlock[i] + cols + 1);
        zmax[i-1] = log(z_temp);
        current_max += zmax[i-1];
        #pragma omp simd
        for (size_t j = 1; j <= cols; ++j) {
            zmBlock[i][j] /= z_temp;
            zeBlock[j] /= z_temp;
            zfBlock[j] /= z_temp;
        }
        
        
        zmBlock[i][0] -= zmax[i-1];
        ze_first[i] -= zmax[i-1];
        zf_first[i] -= current_max;
        if (i < rows) {
            zmBlock[i+1][0] -= current_max;
            ze_first[i+1] -= current_max;
            z_init[0][i-1] = log(zmBlock[i][cols]) + current_max;
            z_init[1][i-1] = log(zeBlock[cols]) + current_max;
            z_init[2][i-1] = log(zfBlock[cols]) + current_max;

        }
        
        // for (size_t j = i; j <= rows; ++j) {
        //     zmBlock[j][0] -= zmax[i-1];
        //     zeBlock[j][0] -= zmax[i-1];
        //     zfBlock[j][0] -= zmax[i-1];
        // }

        //free(zm_exp);
    }
    // Debug(Debug::INFO) << '\n';

    //Calculate the cumulative sum of zmax[1:]
    std::vector<double> rescale(rows);
    // std::partial_sum(zmax + 1, zmax + rows + 1, rescale.begin());
    std::partial_sum(zmax, zmax + rows, rescale.begin());

    //Fixme
    // 
    /*for (size_t i = 0; i < rows; ++i) {
        z_init[0][i] = log(zmBlock[i + 1][memcpy_cols]) + rescale[i];
        z_init[1][i] = log(zeBlock[i + 1][memcpy_cols]) + rescale[i];
        z_init[2][i] = log(zfBlock[i + 1][memcpy_cols]) + rescale[i];
    }*/

    for (size_t i = 0; i < rows; ++i) {
        memcpy(zm[i] + start, zmBlock[i+1]+1, memcpy_cols * sizeof(double));
        //memcpy(ze[i] + start, zeBlock[i+1]+1, memcpy_cols * sizeof(float));
        //memcpy(zf[i] + start, zfBlock[i+1]+1, memcpy_cols * sizeof(float));
    }

    //free(zmBlock);
    //free(zeBlock);
    //free(zfBlock);
    //free(exp_ge_arr);
}

void rescaleBlocks(double **matrix, double **scale, size_t rows, size_t length, size_t blocks, size_t targetLen){
    // Function to rescale the values by the maximum in the log space for each block
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = std::min((b + 1) * length, targetLen);
        // size_t end = (b + 1) * length;
        std::vector<double> cumsum(rows);
        std::partial_sum(scale[b], scale[b] + rows, cumsum.begin());
        // DEBUG:: print cumsum vector for each block
        // std::cout << "block " << b << " cumsum: ";
        // for (size_t i = 0; i < rows; ++i) {
        //     std::cout << cumsum[i] << " ";
        // }
        // std::cout << std::endl;

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = start; j < end; ++j) {
                // matrix[i][j] = log(matrix[i][j]) ;// + cumsum[i-1];
                matrix[i][j] = log(matrix[i][j]) + cumsum[i];
            }
        }
    }
}


double** allocateMemory(size_t rows, size_t cols) {
    // Allocate memory for an array of pointers to rows
    double** array = (double**)malloc(rows * sizeof(double*));
    if (array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Allocate memory for each row
    for (int i = 0; i < rows; i++) {
        array[i] = (double*)malloc(cols * sizeof(double));
        if (array[i] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
    }

    return array;
}
//}



int main() {

    //std::string query = "d1ufaa2"; //d1ufaa2
    //string target = "d2c1ia1"; //d2c1ia1

    size_t rows = 404;
    size_t cols = 196;
    std::cout << "rows: " << rows << " cols: " << cols << std::endl;

    double* scoreForward = new double[rows *cols];
    //float** scoreBackward = allocateMemory(rows, cols);
    double* P =  new double[rows *cols];

    // Load the values into scoreForward and scoreBackward matrices
    std::ifstream forwardFile("./S.txt");
    //std::ifstream backwardFile("/home/lasse/Desktop/Projects/FB_martin/scoreBackward.txt");

    if (forwardFile.is_open()) {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                // Load the values from the files
                forwardFile >> scoreForward[i * cols + j];

                if (forwardFile.fail()) {
                std::cerr << "Error reading from file." << std::endl;
                return 1;
            }
            }
        }
        forwardFile.close();
        //backwardFile.close();
    } else {
        std::cout << "Unable to open score files" << std::endl;
    }

    size_t length = 16;
    size_t blocks = (cols / length) + (cols % length != 0);
    std::cout << blocks << std::endl;

    //align(scoreForward, P, rows, cols, -3.5, -0.3, 10, length, blocks);


    free(scoreForward);
    //free(scoreBackward);
    free(P);

    return 0;
}



