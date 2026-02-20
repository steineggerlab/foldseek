//
// Created by Martin Steinegger on 6/7/21.
//

#ifndef FOLDSEEK_GEMMIWRAPPER_H
#define FOLDSEEK_GEMMIWRAPPER_H
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdlib>
#include "structureto3di.h"

class GemmiWrapper {
public:
    enum class Format {
        Detect = 0,
        Pdb,
        Mmcif,
        Mmjson,
        ChemComp,
        Foldcomp,
        Unknown
    };

    enum class CompressionFormat {
        Detect = 0,
        Gzip = 1,
        Zstd = 2
    };

    GemmiWrapper();
    ~GemmiWrapper() {
        if (fixupBuffer) {
            free(fixupBuffer);
        }
    }

    bool loadFromBuffer(
        const char * buffer,
        size_t bufferSize,
        const std::string& name,
        Format format = Format::Detect,
        CompressionFormat compressionFormat = CompressionFormat::Detect
    );

    bool load(
        const std::string& filename,
        Format format = Format::Detect,
        CompressionFormat compressionFormat = CompressionFormat::Detect
    );

    std::vector<Vec3> ca;
    std::vector<float> ca_bfactor;
    std::vector<Vec3> n;
    std::vector<Vec3> c;
    std::vector<Vec3> cb;
    std::vector<char> ami;
    std::vector<char> seq3di;
    std::vector<std::string> names;
    std::vector<std::string> chainNames;
    std::vector<int> chainStartSerial;
    std::vector<int> chainStartResId;
    std::vector<std::string> chainDescriptions;
    std::vector<unsigned int> modelIndices;
    unsigned int modelCount = 0;
    std::vector<std::pair<size_t, size_t>> chain;
    std::vector<int> taxIds;
    std::string title;

    char* fixupBuffer;
    size_t fixupBufferSize;

private:
    int modelIt;
    int chainIt;

    bool loadFoldcompStructure(std::istream& stream, const std::string& filename);
    void updateStructure(
        void * structure,
        const std::string & filename,
        std::unordered_map<std::string, int>& entity_to_tax_id,
        std::unordered_map<std::string, std::string>& entity_to_description
    );
};

bool GemmiToFoldcomp(
    const GemmiWrapper& gw,
    size_t chainIndex,
    std::string& outBlob,
    int anchorResidueThreshold = 25
);

#endif //FOLDSEEK_GEMMIWRAPPER_H
