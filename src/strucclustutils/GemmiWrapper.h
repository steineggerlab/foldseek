//
// Created by Martin Steinegger on 6/7/21.
//

#ifndef FOLDSEEK_GEMMIWRAPPER_H
#define FOLDSEEK_GEMMIWRAPPER_H
#include <vector>
#include <unordered_map>
#include <string>
#include "structureto3di.h"

class GemmiWrapper {
public:
    GemmiWrapper();

    bool loadFromBuffer(const char * buffer, size_t bufferSize, std::string & name);

    bool load(std::string & filename);

    std::pair<size_t, size_t> nextChain();

    std::vector<Vec3> ca;
    std::vector<float> ca_bfactor;
    std::vector<Vec3> n;
    std::vector<Vec3> c;
    std::vector<Vec3> cb;
    std::vector<char> ami;
    std::vector<std::string> names;
    std::vector<std::string> chainNames;
    std::vector<std::pair<size_t ,size_t>> chain;
    std::string title;
private:
    std::unordered_map<std::string,char> threeAA2oneAA;
    int modelIt;
    int chainIt;

    void updateStructure(void * structure, std::string & filename);
};


#endif //FOLDSEEK_GEMMIWRAPPER_H
