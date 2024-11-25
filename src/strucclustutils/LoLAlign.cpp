#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <vector>
#include "LoLAlign.h"
#include "Fwbw.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "SubstitutionMatrix.h"
#include "Matcher.h"
#include "Util.h"
#include "LocalParameters.h"



int lolalign(int argc, const char **argv, const Command &command) {
    std::cout << "lolAlign" << std::endl;    
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    
    DBReader<unsigned int> qdbrAA(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qdbrAA.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> qdbr3Di((par.db1 + "_ss").c_str(), (par.db1 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qdbr3Di.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> tdbrAA(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tdbrAA.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> tdbr3Di((par.db2 + "_ss").c_str(), (par.db2 + "_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tdbr3Di.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> alnRes (par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    alnRes.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter fwbwAlnWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, LocalParameters::DBTYPE_ALIGNMENT_RES);
    fwbwAlnWriter.open();

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    float aaFactor = (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA) ? 1.4 : 0.0;
    SubstitutionMatrix subMatAA(blosum.c_str(), aaFactor, par.scoreBias);  
    return 1;
}
