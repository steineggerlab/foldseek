#include "ReducedMatrix.h"
#include "Indexer.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"

#include "Debug.h"
#include "AminoAcidLookupTables.h"
#include "DBReader.h"
#include "DBWriter.h"

#include "LocalParameters.h"

#include "kerasify/keras_model.h"
//#include "predict_coding_acc9540_57x32x64.model.h"
//#include "predict_coding_acc9642_57x32x64.model.h"
//#include "predict_coding_acc9598_57x32x64.model.h"
#include "predict_coding_acc9743_57x32x64.model.h"

//#include "predict_coding_acc9260_56x96.model.h"
//#include "predict_coding_acc9623_57x32x64.model.h"

#ifdef OPENMP
#include <omp.h>
#endif

int filternoncoding(int argc, const char **argv, const Command& command)  {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, 2);

    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> seqDb (par.db1.c_str(), par.db1Index.c_str());
    seqDb.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Output file: " << par.db2 << "\n";
    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads));
    dbw.open();

    // Initialize model.
    KerasModel model;
//    model.LoadModel(std::string((const char *)predict_coding_acc9623_57x32x64_model, predict_coding_acc9623_57x32x64_model_len));
    model.LoadModel(std::string((const char *)predict_coding_acc9743_57x32x64_model, predict_coding_acc9743_57x32x64_model_len));

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
//    ReducedMatrix redMat3(subMat.probMatrix, subMat.subMatrixPseudoCounts, 3, subMat.getBitFactor());

    // Create a 1D Tensor on length 20 for input data.
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Sequence seq(par.maxSeqLen, Sequence::AMINO_ACIDS, &subMat,  par.kmerSize, false, false);

#pragma omp for schedule(static)
        for (size_t id = 0; id < seqDb.getSize(); id++) {
            std::vector<float> data;
            char *seqData = seqDb.getData(id);
            unsigned int dbKey = seqDb.getDbKey(id);
            seq.mapSequence(id, dbKey, seqData);

            if(fgets(line,255,typein)==NULL)
                break;
            if(line[11]==' ')
                line_chain='-';
            else
                line_chain=line[11];
            if(chain==line_chain)
            {
                strncpy(ftmp,&line[91],6);
                ftmp[6]='\0';
                kappa=atof(ftmp);
                kappa_tmp=atof(ftmp);
                strncpy(ftmp,&line[97],6);
                ftmp[6]='\0';
                alpha=atof(ftmp);
                alpha_tmp=atof(ftmp);

                if(kappa==180.0)	kappa=179.0;
                if(alpha==180.0)	alpha=179.0;
                kappa=floor(kappa/10);
                alpha=floor((alpha+180)/10);

                if (kappa>18 || alpha>36)
                    continue;
                else
                {
                    output[0]=coding[(int)kappa][(int)alpha];
                }
                if(output[0]=='A')
                    if(kappa_tmp<=114 && alpha_tmp>46)	output[0]='Y';
                strncat(SA_string,output,1);
            }
        }

            dbw.writeData(seqData, seqDb.getSeqLens(id) - 1, dbKey, thread_idx);
        }

        delete[] diAACnt;
//        delete[] pentaAACnt;
    }
//    std::cout << "Filtered: " << static_cast<float>(cnt)/ static_cast<float>(seqDb.getSize()) << std::endl;
    dbw.close(Sequence::AMINO_ACIDS);
    seqDb.close();

    return EXIT_SUCCESS;
}

