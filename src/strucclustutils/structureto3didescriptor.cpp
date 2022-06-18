//
// Created by Martin Steinegger on 11/14/20.
//
#include "Command.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"
#include "FastSort.h"

#include "structureto3di.h"
#include "SubstitutionMatrix.h"
#include "GemmiWrapper.h"
#include "PulchraWrapper.h"

#include <iostream>
#include <dirent.h>

#ifdef OPENMP
#include <omp.h>
#endif

int structureto3didescriptor(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_COMMON);

    std::vector<std::string> filenames(par.filenames);
    std::string outputName = filenames.back();
    filenames.pop_back();
    if(filenames.size() == 1 && FileUtil::directoryExists(filenames.back().c_str())){
        std::vector<std::string> dirs;
        dirs.push_back(filenames.back());
        filenames.pop_back();
        while(dirs.size() != 0){
            std::string dir = dirs.back();
            dirs.pop_back();
            DIR * dpdf = opendir(dir.c_str());
            if (dpdf != NULL) {
                while (dirent * epdf = readdir(dpdf)) {
                    std::string filename(epdf->d_name);
                    if(filename != "." && filename !=".."){
                        if (epdf->d_type == DT_DIR){
                            dirs.push_back(dir+"/"+filename);
                        }else{
                            filenames.push_back(dir+"/"+filename);
                        }
                    }
                }
                closedir(dpdf);
            }
        }
    }
    Debug(Debug::INFO) << "Output file: " << par.db2 << "\n";
    SORT_PARALLEL(filenames.begin(), filenames.end());

    DBWriter vec3di((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    vec3di.open();

    SubstitutionMatrix mat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    Debug::Progress progress(filenames.size());
    std::vector<std::pair<size_t, size_t>> fileIdLookup(filenames.size());
    size_t globalCnt = 0;
    size_t incorrectFiles = 0;
    size_t toShort = 0;

    //===================== single_process ===================//__110710__//
#pragma omp parallel default(none) shared(par, vec3di, mat, filenames, progress, globalCnt, fileIdLookup) reduction(+:incorrectFiles, toShort)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        //recon_related
        StructureTo3Di structureTo3Di;
        PulchraWrapper pulchra;
        GemmiWrapper readStructure;
        std::vector<char> alphabet3di;
        std::vector<float> camol;
        std::string header;
        std::string name;
        std::string result;

#pragma omp for schedule(static)
        for (size_t i = 0; i < filenames.size(); i++) {
            progress.updateProgress();
            // clear memory
            if (readStructure.load(filenames[i]) == false) {
                incorrectFiles++;
                continue;
            }
            size_t id = __sync_fetch_and_add(&globalCnt, readStructure.chain.size());
            fileIdLookup[i].first = id;
            fileIdLookup[i].second = id + readStructure.chain.size();

            for (size_t ch = 0; ch < readStructure.chain.size(); ch++) {
                size_t dbKey = id + ch;
                size_t chainStart = readStructure.chain[ch].first;
                size_t chainEnd = readStructure.chain[ch].second;
                size_t chainLen = chainEnd - chainStart;
                if (chainLen <= 3) {
                    toShort++;
                    continue;
                }

                // Detect if structure is Ca only
                if (std::isnan(readStructure.n[chainStart + 0].x) &&
                    std::isnan(readStructure.n[chainStart + 1].x) &&
                    std::isnan(readStructure.n[chainStart + 2].x) &&
                    std::isnan(readStructure.n[chainStart + 3].x) &&
                    std::isnan(readStructure.c[chainStart + 0].x) &&
                    std::isnan(readStructure.c[chainStart + 1].x) &&
                    std::isnan(readStructure.c[chainStart + 2].x) &&
                    std::isnan(readStructure.c[chainStart + 3].x)) {
                    pulchra.rebuildBackbone(&readStructure.ca[chainStart],
                                            &readStructure.n[chainStart],
                                            &readStructure.c[chainStart],
                                            &readStructure.ami[chainStart],
                                            chainLen);
                }

                header.clear();
                header.append(readStructure.names[ch]);
                if(par.chainNameMode == LocalParameters::CHAIN_MODE_ADD ||
                   (par.chainNameMode == LocalParameters::CHAIN_MODE_AUTO && readStructure.names.size() > 1)){
                    header.push_back('_');
                    header.append(readStructure.chainNames[ch]);
                }
                if(readStructure.title.size() > 0){
                    header.push_back(' ');
                    header.append(readStructure.title);
                }

                char * seq3di = structureTo3Di.structure2states(&readStructure.ca[chainStart],
                                                               &readStructure.n[chainStart],
                                                               &readStructure.c[chainStart],
                                                               &readStructure.cb[chainStart],
                                                               chainLen);
                std::vector<StructureTo3Di::Feature> features = structureTo3Di.getFeatures();
                result.clear();
                result.append(header);
                result.push_back('\t');
                for (size_t j = 0; j < chainLen; j++) {
                    result.push_back(readStructure.ami[chainStart+j]);
                }
                result.push_back('\t');
                for (size_t j = 0; j < chainLen; j++) {
                    result.push_back(mat.num2aa[(size_t)seq3di[j]]);
                }
                result.push_back('\t');
                for (size_t j = 0; j < features.size(); j++) {
                    for(size_t f = 0; f < Alphabet3Di::FEATURE_CNT; f++){
                        result.append(SSTR(features[j].f[f]));
                        result.push_back(',');
                    }
                }
                result[result.size() - 1] = '\n';
                vec3di.writeData((const char*)result.data(), result.size(), dbKey, thread_idx, false);
            }
        }
    }
    vec3di.close(true);
    FileUtil::remove((outputName+".index").c_str());
    return 0;
}
