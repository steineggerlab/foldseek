//
// Created by Martin Steinegger on 11/14/20.
//
#include "Command.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"

#include "structureto3di.h"
#include "mmread.hpp"

#include <iostream>
#include <dirent.h>
#include <mmseqs/src/commons/SubstitutionMatrix.h>

#ifdef OPENMP
#include <omp.h>
#endif


int convert2db(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_COMMON);

    std::vector<std::string> filenames(par.filenames);
    std::string outputName = filenames.back();
    filenames.pop_back();
    if(filenames.size() == 1 && FileUtil::directoryExists(filenames.back().c_str())){
        std::string dir = filenames.back();
        filenames.pop_back();
        DIR * dpdf = opendir(dir.c_str());
        if (dpdf != NULL) {
            while (dirent * epdf = readdir(dpdf)) {
                std::string filename(epdf->d_name);
                if(filename != "." && filename !=".."){
                    filenames.push_back(dir+"/"+filename);
                }
            }
        }
        closedir(dpdf);
    }
    Debug(Debug::INFO) << "Output file: " << par.db2 << "\n";

    DBWriter torsiondbw((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_AMINO_ACIDS);
    torsiondbw.open();
    DBWriter hdbw((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdbw.open();
    DBWriter cadbw((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    cadbw.open();
    DBWriter aadbw((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aadbw.open();
    SubstitutionMatrix mat(par.scoringMatrixFile.aminoacids, 2.0, par.scoreBias);

    size_t incorrectFiles = 0;
    //===================== single_process ===================//__110710__//
#pragma omp parallel default(none) shared(par, torsiondbw, hdbw, cadbw, aadbw, range, filenames) reduction(+:incorrectFiles)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        //recon_related
        StructureTo3Di structureTo3Di;
        std::vector<Vec3> ca;
        std::vector<Vec3> cb;
        std::vector<Vec3> n;
        std::vector<Vec3> c;
        std::vector<char> alphabet3di;
        std::vector<char> ami;
        std::vector<float> camol;
        std::unordered_map<std::string,char> threeAA2oneAA = {
                {"ALA",'A'},  {"ARG",'R'},  {"ASN",'N'}, {"ASP",'D'},
                {"CYS",'C'},  {"GLN",'Q'},  {"GLU",'E'}, {"GLY",'G'},
                {"HIS",'H'},  {"ILE",'I'},  {"LEU",'L'}, {"LYS",'K'},
                {"MET",'M'},  {"PHE",'F'},  {"PRO",'P'}, {"SER",'S'},
                {"THR",'T'},  {"TRP",'W'},  {"TYR",'Y'}, {"VAL",'V'}
        };
        std::string name;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < filenames.size(); i++) {

            // clear memory
            gemmi::Structure st = gemmi::read_structure_file(filenames[i], gemmi::CoorFormat::Pdb);
            if(st.models.size() == 0 || st.models[0].chains.size() == 0 ||
               st.models[0].chains[0].residues.size()  == 0 ){
                incorrectFiles++;
                continue;
            }
            for (gemmi::Model& model : st.models){
                for (gemmi::Chain& chain : model.chains){
                    name.clear();
                    ca.clear(); c.clear(); cb.clear(); n.clear(); ami.clear(); alphabet3di.clear();
                    name = FileUtil::baseName(filenames[i]);
                    name.push_back('_');
                    name.append(chain.name);
                    for (gemmi::Residue& res : chain.residues){
                        if (res.het_flag != 'A')
                            continue;
                        Vec3 ca_atom = {NAN, NAN, NAN};
                        Vec3 cb_atom = {NAN, NAN, NAN};
                        Vec3 n_atom = {NAN, NAN, NAN};
                        Vec3 c_atom = {NAN, NAN, NAN};

                        for (gemmi::Atom& atom : res.atoms){
                            if (atom.name == "CA"){
                                ca_atom.x = atom.pos.x;
                                ca_atom.y = atom.pos.y;
                                ca_atom.z = atom.pos.z;
                            }
                            else if (atom.name == "CB"){
                                cb_atom.x = atom.pos.x;
                                cb_atom.y = atom.pos.y;
                                cb_atom.z = atom.pos.z;
                            }
                            else if (atom.name == "N"){
                                n_atom.x = atom.pos.x;
                                n_atom.y = atom.pos.y;
                                n_atom.z = atom.pos.z;
                            }
                            else if (atom.name == "C"){
                                c_atom.x = atom.pos.x;
                                c_atom.y = atom.pos.y;
                                c_atom.z = atom.pos.z;
                            }
                        }
                        ca.push_back(ca_atom);
                        cb.push_back(cb_atom);
                        n.push_back(n_atom);
                        c.push_back(c_atom);
                        ami.push_back(threeAA2oneAA[res.name]);
                    }
                    char * states = structureTo3Di.structure2states(ca, n, c, cb);
                    for(size_t pos = 0; pos < ca.size(); pos++){
                        alphabet3di.push_back(mat.num2aa[static_cast<int>(states[pos])]);
                    }
                    alphabet3di.push_back('\n');

                    torsiondbw.writeData(alphabet3di.data(), alphabet3di.size(), i, thread_idx);
                    ami.push_back('\n');
                    aadbw.writeData(ami.data(), ami.size(), i, thread_idx);
                    name.push_back('\n');
                    hdbw.writeData(name.c_str(), name.size(), i, thread_idx);
                    name.clear();
                    for(size_t res = 0; res < ca.size(); res++){
                        camol.push_back(ca[res].x);
                    }
                    for(size_t res = 0; res < ca.size(); res++) {
                        camol.push_back(ca[res].y);
                    }
                    for(size_t res = 0; res < ca.size(); res++) {
                        camol.push_back(ca[res].z);
                    }
                    cadbw.writeData((const char*)camol.data(), camol.size() * sizeof(float), i, thread_idx);

                    break;
                }
                break;
            }


        }
    }


    torsiondbw.close(true);
    hdbw.close(true);
    cadbw.close(true);
    aadbw.close(true);
    //if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
    DBReader<unsigned int>::softlinkDb(outputName, outputName+"_ca", DBFiles::HEADERS);
    DBReader<unsigned int>::softlinkDb(outputName, outputName+"_ss", DBFiles::HEADERS);
    //}
    Debug(Debug::INFO) << incorrectFiles << " out of " << filenames.size() << " entries are incorrect.\n";
    return EXIT_SUCCESS;
}
