//
// Created by Martin Steinegger on 11/14/20.
//
#include "Command.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"
#include "XYZ.h"
#include "Confo_Lett.h"
#include "Mol_File.h"
#include <iostream>
#include <dirent.h>

#ifdef OPENMP
#include <omp.h>
#endif

using namespace std;
//------------- Get_PDB_File_Len -----------//
int Get_PDB_File_Len(string &pdbfile, string &name) //-> only suitable for pdb_BC100 pdb_file
{
    //--- list for mapping ---//
    map<string, int> ws_mapping;
    map<string, int>::iterator iter;
    //read
    ifstream fin;
    string buf,temp;
    fin.open(pdbfile.c_str(), ios::in);
    if(fin.fail()!=0)
    {
        Debug(Debug::ERROR) << "File "<< pdbfile <<"not found!!\n";
        EXIT(EXIT_FAILURE);
    }
    //process
    int len;
    int count=0;
    const char * words[512];
    bool headerFound = false;
    for(;;)
    {
        //HEADER    ELECTRON TRANSFER(CUPROPROTEIN)         28-JUN-88   1PAZ
        if(!getline(fin,buf,'\n'))break;
        len=(int)buf.length();
        if(len<3)continue;
        //check ATOM
        if(len<4)continue;
        temp=buf.substr(0,6);
        if(temp =="HEADER" ) {
            headerFound = true;
            buf.erase(std::find_if(buf.rbegin(), buf.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), buf.end());
            std::size_t len = Util::getWordsOfLine(buf.c_str(), words, 512);
            size_t keyLen = words[len]-words[len-1];
            name.assign(words[len-1], keyLen);
            continue;
        }


        temp=buf.substr(0,4);
        if(temp!="ATOM" && temp!="HETA")continue;
        //check CA
        temp=buf.substr(13,2);
        if(temp!="CA") continue;
        //record name
        temp=buf.substr(21,6);
        iter = ws_mapping.find(temp);
        if(iter != ws_mapping.end())continue;
        count++;
        ws_mapping.insert(map < string, int >::value_type(temp, count));
    }
    if(headerFound == false){
        name.assign(FileUtil::baseName(pdbfile));
    }
    //return
    return count;
}


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

    DBWriter torsiondbw((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_TOSION_SEQUENCE);
    torsiondbw.open();
    DBWriter hdbw((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdbw.open();
    DBWriter cadbw((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    cadbw.open();
    DBWriter aadbw((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aadbw.open();

    const int INPUT_MODE=1; //main (default:1)
    const int INPUT_TYPE=1; //main (default:1)
    const int INPUT_GLYS=1; //vice (default:1)
    const int WARN_OUT=1;   //vice (default:1)

    string range="_";
    size_t incorrectFiles = 0;
    //===================== single_process ===================//__110710__//
#pragma omp parallel default(none) shared(par, torsiondbw, hdbw, cadbw, aadbw, range, filenames) reduction(+:incorrectFiles)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        int allocSize = static_cast<int>(par.maxSeqLen);
        //recon_related
        XYZ *mol = new XYZ[allocSize]; //CA
        float * camol = new float[allocSize*3]; //CA
        char *cle = new char[allocSize + 1];
        char *ami = new char[allocSize + 1];
        PDB_Residue *pdb = new PDB_Residue[allocSize];
        std::string name;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < filenames.size(); i++) {
            //get length
            int totalLen = Get_PDB_File_Len(filenames[i], name);
            if (totalLen <= 0) {
                incorrectFiles++;
                continue;
            }

            if(totalLen > allocSize){
                //delete
                delete[] pdb;
                delete[] mol;
                delete[] ami;
                delete[] cle;
                delete[] camol;
                //create
                pdb = new PDB_Residue[totalLen];
                mol = new XYZ[totalLen];
                camol = new float[totalLen * 3]; //CA
                ami = new char[totalLen + 1];
                cle = new char[totalLen + 1];
                allocSize = totalLen;
            }

            //class
            Mol_File mol_input;
            Confo_Lett confo_lett;
            mol_input.MODRES = 1;
            //init
            int moln;
            string TER = "TER                                                                             ";
            //open
            mol_input.CbBACK = 1;
            mol_input.CaONLY = 0; //consider Non-CA atoms !!//__110408__//
            //macro
            mol_input.OUTP_MODE = INPUT_TYPE;
            mol_input.PROC_MODE = INPUT_MODE;
            mol_input.GLYC_MODE = INPUT_GLYS;
            //memory limit
            mol_input.WARNING_out = WARN_OUT;
            mol_input.MEMORY_LIMIT = totalLen;

            //process
            {
                //pre_process
                int ret_val = mol_input.XYZ_Input(filenames[i], range, 0, moln, 0, 0, 0, 0, 0);

                if (ret_val <= 0) {
                    if (ret_val != -12345){
                        incorrectFiles++;
                        continue;
                    }
                }
                if (ret_val != 1){
                    incorrectFiles++;
                    continue;
                }
                //check memory
                if (ret_val == -12345) {
                    //add memory
                    if (moln > allocSize) {
                        totalLen = moln;
                        allocSize = totalLen;
                        //delete
                        delete[] pdb;
                        delete[] mol;
                        delete[] ami;
                        delete[] cle;
                        delete[] camol;
                        //create
                        pdb = new PDB_Residue[totalLen];
                        mol = new XYZ[totalLen];
                        ami = new char[totalLen + 1];
                        cle = new char[totalLen + 1];
                        camol = new float[totalLen * 3]; //CA
                        //memory limit
                        mol_input.MEMORY_LIMIT = totalLen;
                    }
                }
                //reload
                int PRE_LOAD_ = mol_input.PRE_LOAD;
                int WARNING_out_ = mol_input.WARNING_out;
                mol_input.PRE_LOAD = 1;
                mol_input.WARNING_out = 0;
                mol_input.XYZ_Input(filenames[i], range, 0, totalLen, mol, ami, 0, 0, pdb);
                mol_input.PRE_LOAD = PRE_LOAD_;
                mol_input.WARNING_out = WARNING_out_;

                //check if PDB files was correct
                {
                    int correct = 1; //default:OK
                    for (int res = 0; res < totalLen; res++) {
                        int iret = pdb[res].PDB_residue_backbone_check(4);
                        if (iret != 1) {
                            correct = 0;
                            break;
                        }
                        iret = pdb[res].PDB_residue_CB_check();
                        if (iret != 1) {
                            correct = 0;
                            break;
                        }
                    }
                    if(correct == 0){
                        incorrectFiles++;
                        continue;
                    }

                    confo_lett.btb_ori(0, 0, 0, totalLen, mol, cle);
                    cle[totalLen]='\n';
                    torsiondbw.writeData(cle, totalLen + 1, i, thread_idx);
                    ami[totalLen]='\n';
                    aadbw.writeData(ami, totalLen + 1, i, thread_idx);
                    name.push_back('\n');
                    hdbw.writeData(name.c_str(), name.size(), i, thread_idx);
                    name.clear();
                    for(int res = 0; res < totalLen; res++){
                        camol[res*3 + 0] = mol[res].X;
                        camol[res*3 + 1] = mol[res].Y;
                        camol[res*3 + 2] = mol[res].Z;
                    }
                    cadbw.writeData((const char*)camol, totalLen * 3 * sizeof(float), i, thread_idx);
                }
            }
        }
        //delete
        delete[] pdb;
        delete[] mol;
        delete[] ami;
        delete[] cle;
        delete[] camol;
    }
    torsiondbw.close(true);
    hdbw.close(true);
    cadbw.close(true);
    aadbw.close(true);
    //if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
    DBReader<unsigned int>::softlinkDb(outputName, outputName+"_ca", DBFiles::HEADERS);
    DBReader<unsigned int>::softlinkDb(outputName, outputName+"_tortion", DBFiles::HEADERS);
    //}
    Debug(Debug::INFO) << incorrectFiles << " out of " << filenames.size() << " entries are incorrect.\n";
    return EXIT_SUCCESS;
}
