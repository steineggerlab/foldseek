
#include <string>
#include <vector>
#include <tmalign/TMalign.h>
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "LocalParameters.h"
#include "Matcher.h"

#ifdef OPENMP
#include <omp.h>
#endif


int tmalign(int argc, const char **argv, const Command& command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, 4);

    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tdbr = NULL;

    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.noPreload == false) {
            tdbr->readMmapedDataInMemory();
        }
    }


    Debug(Debug::INFO) << "Result database: " << par.db3 << "\n";
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);


    Debug(Debug::INFO) << "Output file: " << par.db4 << "\n";
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads));
    dbw.open();



    bool i_opt = false; // flag for -i, with user given initial alignment
    bool I_opt = false; // flag for -I, stick to user given alignment
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    bool fast_opt = false; // flags for -fast, fTM-align algorithm
    double Lnorm_ass = 0.0;
    double  d0_scale = 0.0;


    // parsing code
    string atom_opt = "auto";// use C alpha atom for protein and C3' for RNA
    string suffix_opt = "";    // set -suffix to empty
    string dir_opt = "";    // set -dir to empty
    string dir1_opt = "";    // set -dir1 to empty
    string dir2_opt = "";    // set -dir2 to empty
    int ter_opt = 1;     // TER, END, or different chainID
    int split_opt = 2;     // split chain
    int infmt1_opt = 0;     // PDB format for chain_1
    int infmt2_opt = 0;
    int byresi_opt = 0;     // set -byresi to 0
#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string seqM, seqxA, seqyA;// for output alignment
        vector<vector<string> > PDB_lines1; // text of chain1
        vector<vector<string> > PDB_lines2; // text of chain2
        // ya[0...ylen-1][0..2], in general,
        // ya is regarded as native structure
        // --> superpose xa onto ya
        std::vector<std::string> resi_vec1;  // residue index for chain1
        std::vector<std::string> resi_vec2;  // residue index for chain2
        std::vector<std::string> chainID_list1;
        std::vector<std::string> chainID_list2;
        std::vector<std::string> sequence; // get value from alignment file
        char *seqx, *seqy;       // for the protein sequence
        int *secx, *secy;       // for the secondary structure
        int *xresno, *yresno;   // residue number for fragment gapless threading
        double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
/* declare variable specific to this pair of TMalign */
        double t0[3], u0[3][3];
        double TM1, TM2;
        double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
        double d0_0, TM_0;
        double d0A, d0B, d0u, d0a;
        double d0_out = 5.0;
        double rmsd0 = 0.0;
        int L_ali;                // Aligned length in standard_TMscore
        double Liden = 0;
        double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
        int n_ali = 0;
        int n_ali8 = 0;
        char buffer[1024+32768];
        std::string resultBuffer;
#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            Debug::printProgress(id);

            char *data = resultReader.getData(id);
            size_t queryKey = resultReader.getDbKey(id);
            unsigned int queryId = qdbr.getId(queryKey);
            char *querySeq = qdbr.getData(queryId);
            unsigned int chain_i = 0;

            std::istringstream xpdb(querySeq);


            int numChains = get_PDB_lines(xpdb, PDB_lines1, chainID_list1,
                                          resi_vec1, byresi_opt, ter_opt, infmt1_opt, atom_opt, split_opt);

            int xlen = PDB_lines1[chain_i].size();
            if (!xlen) {
                Debug(Debug::ERROR) << "Warning! Cannot parse file: " << queryKey
                                    << ". Chain length 0.\n";
                continue;
            } else if (xlen <= 5) {
                Debug(Debug::ERROR) << "Sequence is too short <=5! Model:  " << queryKey << "\n";
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new int[xlen];
            xresno = new int[xlen];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, xresno);
            make_sec(xa, xlen, secx); // secondary structure assignment
            std::vector<hit_t> results = QueryMatcher::parsePrefilterHits(data);
            for (size_t entryIdx = 0; entryIdx < results.size(); entryIdx++) {
                unsigned int targetId = tdbr->getId(results[entryIdx].seqId);
                const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;
                if(isIdentity == true){
                    std::string backtrace = "";
                    Matcher::result_t result(results[entryIdx].seqId, 0 , 1.0, 1.0, 1.0, TM1, std::max(xlen,xlen), 0, xlen-1, xlen, 0, xlen-1, xlen, backtrace);
                    size_t len = Matcher::resultToBuffer(buffer, result, true, false);
                    resultBuffer.append(buffer, len);
                    continue;
                }
                char * targetSeq = tdbr->getData(targetId);
                unsigned int chain_j = 0;
                std::istringstream ypdb(targetSeq);

                int targetNumChains = get_PDB_lines(ypdb, PDB_lines2, chainID_list2,
                                                    resi_vec2, byresi_opt, ter_opt, infmt2_opt, atom_opt, split_opt);

                int ylen = PDB_lines2[chain_j].size();

                if(Util::canBeCovered(par.covThr, par.covMode, xlen, ylen)==false){
                    continue;
                }

                if (!ylen) {
                    Debug(Debug::ERROR) << "Warning! Cannot parse file: " << results[entryIdx].seqId
                                        << ". Chain length 0.\n";
                    continue;
                } else if (ylen <= 5) {
                    Debug(Debug::ERROR) << "Sequence is too short <=5! Model: " << results[entryIdx].seqId << "\n";
                    continue;
                }
                NewArray(&ya, ylen, 3);
                seqy = new char[ylen + 1];
                yresno = new int[ylen];
                secy = new int[ylen];

                ylen = read_PDB(PDB_lines2[chain_j], ya, seqy, yresno);

                make_sec(ya, ylen, secy); // secondary structure assignment
                /* entry function for structure alignment */
                TMalign_main(
                        xa, ya, xresno, yresno, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, I_opt, a_opt, u_opt, d_opt, fast_opt);
                double seqId = Liden/(static_cast<double>(n_ali8));
                int rmsdScore = static_cast<int>(rmsd0*1000.0);
                std::string backtrace = "";
                Matcher::result_t result(results[entryIdx].seqId, rmsdScore , 1.0, 1.0, seqId, TM1, std::max(xlen,ylen), 0, xlen-1, xlen, 0, ylen-1, ylen, backtrace);


                bool hasCov = Util::hasCoverage(par.covThr, par.covMode, 1.0, 1.0);
                bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                bool hasTMscore = (TM1 >= par.tmScoreThr);
                if(hasCov && hasSeqId  && hasTMscore){
                    size_t len = Matcher::resultToBuffer(buffer, result, true, false);
                    resultBuffer.append(buffer, len);
                }

                DeleteArray(&ya, ylen);

                chainID_list2.clear();
                PDB_lines2.clear();
                resi_vec2.clear();
                sequence.clear();
                delete[] seqy;
                delete[] secy;
                delete[] yresno;
            }
            dbw.writeData(resultBuffer.c_str(), resultBuffer.size(), queryKey, thread_idx);
            resultBuffer.clear();
            chainID_list1.clear();
            PDB_lines1.clear();
            resi_vec1.clear();
            DeleteArray(&xa, xlen);
            delete[] seqx;
            delete[] secx;
            delete[] xresno;
        }
    }

    dbw.close();
    resultReader.close();
    qdbr.close();
    if(sameDB == false){
        tdbr->close();
    }
    return EXIT_SUCCESS;
}


