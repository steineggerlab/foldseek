#include <dssp/structure.h>
#include "StructSubstitutionMatrix.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"

#include "LocalParameters.h"


#ifdef OPENMP
#include <omp.h>
#endif

int convert2statealphabet(int argc, const char **argv, const Command& command)  {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_COMMON);


    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> dsspDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    dsspDb.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Output file: " << par.db2 << "\n";
    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    dbw.open();


    Debug(Debug::INFO) << "Output file: " << par.db2 << "\n";
    DBWriter dbwh((par.db2+"_h").c_str(), (par.db2+"_h.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    dbwh.open();

    char coding[18][36];
    std::string SA_coding_rule =  "EEFFFFZZZZZZZZZZZZZZZZZZZZZZZZZZNNHHEEEFFFKKZZZZZZZZZZZZZZZZZZZZZUUUUNHHHHHKFF"
                                  "KKKKKZZZZZZZZZZZZZZZZZUUUUHHHHHNNNKKKKKKKKKZZZZZZZZZZZZUZUUUUUHHHHNPTTTNNKKKKK"
                                  "KKUKUUUUZUUUUUUUUUUUUNNNPPPPTTTTNNNTUUUUUUUUUUUUUUUUUUUUURRPPPPPPTTTTTTTTTVVVV"
                                  "VUMUMMMMUUURRRRRPPPPPPPPTTTTTTVVVVVVVMMMMMMQQQQQQRRQPPPPSSSSTTTTTTTVVVVVVVVMGM"
                                  "GGQQQQQQQQRQQSSSSSSSSSSVVVVVVVVVVMGGGGQQQQQQQQZZZZZSSSSSSSSSSWWWVVVVMDDBDLQQQQ"
                                  "QRZZZZZZZZZSSSSSSSSWWWWWWWLDACIILLQRRZZZZZZZZZZZSSSSSSWWWWWWLLLLLIILLLRZZZZZZZ"
                                  "ZZZZZZZZSSSSSWWLWLLLLLIILZZZZZZZZZZZZZZZZZZZZZZZZZZZZZLZZZZZZZZZZZZZZZZZZZZZZZ"
                                  "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
                                  "ZZZZZZZZZZZZZZZZZZZZZZZZ" ;
    char output[1];

    strcpy(output,"");

    int k,i,j;
    k=0;
    for(i=0;i<18;i++)
    {
        for(j=0;j<36;j++)
        {
            coding[i][j]=SA_coding_rule[k];
            k++;
        }
    }


    int total = 0;
#pragma omp parallel shared(total)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string headerString;
        std::string outputString;

#pragma omp for schedule(static)
        for (size_t id = 0; id < dsspDb.getSize(); id++) {
            MProtein a;
            char *data = dsspDb.getData(id, thread_idx);
            std::istringstream pdb(data);
            a.ReadPDB(pdb);
            a.CalculateSecondaryStructure();
            unsigned int dbKey = dsspDb.getDbKey(id);
            for(size_t chainIdx = 0; chainIdx < a.GetChains().size(); chainIdx++){
                MChain * chain = a.GetChains()[chainIdx];
                for(size_t resIdx = 0; resIdx < chain->GetResidues().size(); resIdx++) {
                    MResidue * res =  chain->GetResidues()[resIdx];
                    float kappa = res->Kappa();
                    float kappa_tmp = kappa;
                    float alpha = res->Alpha().first;
                    float alpha_tmp = alpha;

                    if (kappa == 180.0) kappa = 179.0;
                    if (alpha == 180.0) alpha = 179.0;
                    kappa = floor(kappa / 10);
                    alpha = floor((alpha + 180) / 10);

                    if (kappa > 18 || alpha > 36) {
                        continue;
                    } else {
                        output[0] = coding[(int) kappa][(int) alpha];
                    }
                    if (output[0] == 'A') {
                        if (kappa_tmp <= 114 && alpha_tmp > 46) {
                            output[0] = 'Y';
                        }
                    }
                    outputString.push_back(output[0]);
                }
                size_t id = __sync_fetch_and_add(&total, 1);

                outputString.push_back('\n');
                dbw.writeData(outputString.c_str(), outputString.size(), id, thread_idx);

                headerString.append(SSTR(dbKey));
                headerString.push_back(':');
                headerString.append(SSTR(chainIdx));
                headerString.push_back('\n');
                dbwh.writeData(headerString.c_str(), headerString.size(), id, thread_idx);

                outputString.clear();
                headerString.clear();
            }

        }
    }
    dbw.close();
    dbwh.close();
    dsspDb.close();
    return EXIT_SUCCESS;
}

