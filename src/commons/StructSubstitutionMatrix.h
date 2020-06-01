#ifndef STRUCTSUBSTITUTION_MATRIX_H
#define STRUCTSUBSTITUTION_MATRIX_H

#include <cstddef>
#include "BaseMatrix.h"
#include "ProfileStates.h"

class StructSubstitutionMatrix: public BaseMatrix {

    public:
        StructSubstitutionMatrix(const char *scoringMatrixFileName_, float bitFactor, float scoreBias);

        virtual ~StructSubstitutionMatrix();

        virtual float getBitFactor() {return bitFactor; }
    
        virtual double getBackgroundProb(size_t aa_index) { return pBack[aa_index]; }

        static void calcLocalAaBiasCorrection(const BaseMatrix *m ,const int *int_sequence, const int N, float *compositionBias);
        static void calcProfileProfileLocalAaBiasCorrection(short *profileScores,
                                                const size_t profileAASize,
                                                const int N,
                                                size_t alphabetSize);
        static void calcProfileProfileLocalAaBiasCorrectionAln(int8_t *profileScores,
                                                             int N,
                                                             size_t alphabetSize,
                                                             BaseMatrix *subMat);
        static void calcGlobalAaBiasCorrection(const BaseMatrix * m,
                                               short *profileScores,
                                               float *pNullBuffer,
                                               const size_t profileAASize,
                                               const int N);
        bool estimateLambdaAndBackground(const double ** mat, int alphabetSize, double * pBack, double & lambda);


        void setupLetterMapping();


        struct FastMatrix{
            const char ** matrix;
            const char * matrixData;
            const size_t asciiStart;
            FastMatrix(const char ** matrix, const char * matrixData, const size_t asciiStart):
                    matrix(matrix), matrixData(matrixData), asciiStart(asciiStart)
            {}
        };

        // build matrix from ~ (=0) to ~(=122)
        static FastMatrix createAsciiSubMat(BaseMatrix & submat){
            const size_t asciiStart = 0;
            const size_t asciiEnd = 'z'+1;
            const size_t range = asciiEnd-asciiStart;
            char ** matrix = new char *[range];
            char * matrixData = new char[range*range];
            for(size_t i = 0; i < range; i++) {
                matrix[i] = matrixData+(i*range);
                int curr_i = submat.aa2int[asciiStart+i];
                for (size_t j = 0; j < range; j++) {
                    int curr_j = submat.aa2int[asciiStart+j];
                    matrix[i][j] = static_cast<char>(submat.subMatrix[curr_i][curr_j]);
                }
            }
            return FastMatrix((const char**) matrix,
                              (const char*) matrixData,
                              asciiStart);
        }

        int alphabetSize;

private:

        const char* scoringMatrixFileName;

        int parseAlphabet(char * word, char * int2aa, int * aa2int);

        int readProbMatrix(const std::string &matrixData);

        float bitFactor;
};

#endif
