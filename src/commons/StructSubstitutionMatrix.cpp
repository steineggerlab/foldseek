#include "StructSubstitutionMatrix.h"
#include "Util.h"
#include "Debug.h"
#include "MathUtil.h"
#include "lambda_calculator.h"
#include <cstring>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <climits>
#include "BLAST3D.out.h"



StructSubstitutionMatrix::StructSubstitutionMatrix(const char *scoringMatrixFileName_,
                                                   float bitFactor, float scoreBias = -0.2) :
        scoringMatrixFileName(scoringMatrixFileName_), alphabetSize(24)
{
    setupLetterMapping();

    std::string submat((const char *)BLAST3D_out,BLAST3D_out_len);
    matrixName = "blast3d.out";
    int alphabetSize = readProbMatrix(submat);
    if (alphabetSize < this->alphabetSize - 1) {
        this->alphabetSize = alphabetSize;
    }

    //print(probMatrix, int2aa, alphabetSize);
    generateSubMatrix(this->probMatrix, this->subMatrixPseudoCounts,
                      this->subMatrix, this->subMatrix2Bit,
                      this->alphabetSize, true, bitFactor, scoreBias);
    this->bitFactor = bitFactor;
}


bool StructSubstitutionMatrix::estimateLambdaAndBackground(const double **scoreMatrix,
                                                           int alphabetSize, double *pBack, double &lambda) {
    // We need to pass the parameters as 1-based pointers, hence the +1s
    // and -1s.

    std::vector<double> cells(alphabetSize * (alphabetSize + 1));
    std::vector<const double *> pointers(alphabetSize + 1);

    for (int i = 0; i < alphabetSize; ++i) {
        pointers[i + 1] = &cells[i * alphabetSize];
        for (int j = 0; j < alphabetSize; ++j) {
            cells[i * alphabetSize + j + 1] = scoreMatrix[i][j];
        }
    }

    std::vector<double> letterProbs1(alphabetSize, 0);
    std::vector<double> letterProbs2(alphabetSize, 0);


    lambda = calculate_lambda(&pointers[0], alphabetSize,
                              &letterProbs1[0] - 1,
                              &letterProbs2[0] - 1);

    for (int i = 0; i < alphabetSize; i++) {
        pBack[i] = letterProbs1[i];
    }

    if (lambda < 0)
        return false; //bad
    else
        return true; //good
}



StructSubstitutionMatrix::~StructSubstitutionMatrix() {
}

void StructSubstitutionMatrix::setupLetterMapping(){
    for(int letter = 0; letter < UCHAR_MAX; letter++) {
        char upperLetter = toupper(static_cast<char>(letter));
        switch (upperLetter){

            case 'A':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'G':
            case 'H':
            case 'I':
            case 'K':
            case 'L':
            case 'M':
            case 'N':
            case 'P':
            case 'Q':
            case 'R':
            case 'S':
            case 'T':
            case 'V':
            case 'W':
            case 'Y':
            case 'X':
            case 'Z':
            case '[':
            case '\\':
                this->aa2int[static_cast<int>(letter)] = this->aa2int[static_cast<int>(upperLetter)];
                break;
            default:
                this->aa2int[static_cast<int>(letter)] = this->aa2int[(int) 'X'];
                break;
        }
    }
}

int StructSubstitutionMatrix::parseAlphabet(char *word, char *int2aa, int *aa2int) {
    char *charReader = word;
    int minAAInt = INT_MAX;
    // find amino acid with minimal int value
    while (isgraph(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2int[static_cast<int>(aa)];
        minAAInt = std::min(minAAInt, intAA);
        charReader++;
    }
    char minAAChar = int2aa[minAAInt];
    // do alphbet reduction
    charReader = word;
    while (isgraph(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2int[static_cast<int>(aa)];
        aa2int[static_cast<int>(aa)] = minAAInt;
        int2aa[intAA] = minAAChar;
        charReader++;
    }
    return minAAInt;
}

int StructSubstitutionMatrix::readProbMatrix(const std::string &matrixData) {
    std::stringstream in(matrixData);
    std::string line;
    bool probMatrixStart = false;

    char *words[256];
    int column_aa[32];
    int column_aa_sorted[32];
    int alphabetSize = 0;
    bool hasLambda = false;
    bool hasBackground = false;
    while (in.good()) {
        getline(in, line);
        size_t wordCnt = Util::getWordsOfLine((char *) line.c_str(), words, 256);
        // skip comments
        if (line[0] == '#') {
            if (line.find("# Background (precomputed optional):") == 0) {
                for (size_t i = 4; i < wordCnt; i++) {
                    double f = strtod(words[i], NULL);
                    pBack[i-4] = f;
                }
                hasBackground = true;
            }
            if (line.find("# Lambda     (precomputed optional):") == 0) {
                double f = strtod(words[4], NULL);
                lambda = f;
                hasLambda = true;
            }
            continue;
        }
        if (wordCnt > 1 && probMatrixStart == false) {
            for (size_t i = 0; i < wordCnt; i++) {
                if (isgraph(words[i][0]) == false) {
                    Debug(Debug::ERROR) << "Probability matrix must start with alphabet header.\n";
                    EXIT(EXIT_FAILURE);
                }
                column_aa[i] = parseAlphabet(words[i], int2aa, aa2int);
            }
            alphabetSize = wordCnt;
            if (alphabetSize > 32) {
                Debug(Debug::ERROR) << "Only alphabets with up to 32 letters are allowed.\n";
                EXIT(EXIT_FAILURE);
            }
            memcpy(column_aa_sorted, column_aa, sizeof(int) * alphabetSize);
            std::sort(column_aa_sorted, column_aa_sorted + alphabetSize);
            int column_old_aa[32];
            memcpy(column_old_aa, column_aa, sizeof(int) * alphabetSize);
            std::map<int, int> mapping;
            for (int i = 0; i < alphabetSize; i++) {
                for (int j = 0; j < alphabetSize; j++) {
                    if (column_aa_sorted[i] == column_aa[j]) {
                        const char repAA = int2aa[column_aa[j]];
                        for (size_t z = 'A'; z < 'Z'; z++) {
                            aa2int[z] = (aa2int[z] == column_aa_sorted[i]) ? i : aa2int[z];
                        }
                        int2aa[i] = repAA;
                        column_aa[j] = i;
                    }
                }
            }
            probMatrixStart = true;
            continue;
        }
        if (wordCnt > 1 && probMatrixStart == true) {
            if (isgraph(words[0][0]) == false) {
                Debug(Debug::ERROR) << "First element in probability line must be an alphabet letter.\n";
                EXIT(EXIT_FAILURE);
            }
            int aa = parseAlphabet(words[0], int2aa, aa2int);
            for (int i = 0; i < alphabetSize; i++) {
                double f = strtod(words[i + 1], NULL);
                probMatrix[aa][column_aa[i]] = f; // divided by 2 because we scale bit/2 ot bit
            }
        }
    }
    bool containsX = false;
    bool xIsPositive = false;
    for (int i = 0; i < alphabetSize; i++) {
        if (column_aa[i] == aa2int[(int)'X']) {
            containsX = true;
            for (int j = 0; j < alphabetSize; j++) {
                int xIndex = aa2int[(int)'X'];
                if ((probMatrix[xIndex][j] > 0) || (probMatrix[j][xIndex] > 0)) {
                    xIsPositive = true;
                    break;
                }
            }
            break;
        }
    }
    if (containsX == false) {
        Debug(Debug::ERROR) << "Please add X to your substitution matrix.\n";
        EXIT(EXIT_FAILURE);
    }

    if(hasLambda == false || hasBackground == false){
        if (estimateLambdaAndBackground(const_cast<const double **>(probMatrix), alphabetSize - ((xIsPositive) ? 0 : 1),
                                        pBack, lambda) == false) {
            Debug(Debug::ERROR) << "Computing inverse of substitution matrix failed\n";
            EXIT(EXIT_FAILURE);
        }
        pBack[aa2int[(int)'\\']]=ANY_BACK;
    }
    if(xIsPositive == false){
        for (int i = 0; i < alphabetSize - 1; i++) {
            pBack[i] = pBack[i] * (1.0 - pBack[aa2int[(int)'\\']]);
        }
    }
    // Reconstruct Probability Sab=(Pab/Pa*Pb) -> Pab = exp(Sab) * Pa * Pb
    for (int i = 0; i < alphabetSize; i++) {
        //smat[i] = smatData+((subMat.alphabetSize-1)*i);
        for (int j = 0; j < alphabetSize; j++) {
            probMatrix[i][j] = std::exp(lambda * probMatrix[i][j]) * pBack[i] * pBack[j];
        }
    }

    return alphabetSize;
}



