//
// Created by Martin Steinegger on 6/7/21.
//

#ifndef FOLDSEEK_GEMMIWRAPPER_H
#define FOLDSEEK_GEMMIWRAPPER_H
#include <vector>
#include <iosfwd>
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
        bool saveResIndex = false,
        Format format = Format::Detect,
        CompressionFormat compressionFormat = CompressionFormat::Detect
    );

    bool load(
        const std::string& filename,
        bool saveResIndex = false,
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
    std::vector<unsigned int> resIds;
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
        std::unordered_map<std::string, std::string>& entity_to_description,
        bool saveResIndex
    );
};

bool GemmiToFoldcomp(
    const GemmiWrapper& gw,
    size_t chainIndex,
    std::string& outBlob,
    int anchorResidueThreshold = 25
);

// Writes structures with gemmi. It lives in this library because gemmi needs
// exceptions, which only the gemmiwrapper target is compiled with, so the
// interface here stays free of gemmi types.
//
// Usage: reset(), then addModel()/addChain()/addCalpha() in that order, then
// writeCif(). The writer can be reused after writeCif().
class GemmiWriter {
public:
    GemmiWriter();
    ~GemmiWriter();

    // drops everything accumulated so far and starts a new entry
    void reset(const std::string& entryName, const std::string& title);
    // starts a new model, modelNumber becomes _atom_site.pdbx_PDB_model_num.
    // mmCIF has no per model title item, a non empty modelTitle is kept in a
    // _foldseek_model loop instead.
    void addModel(int modelNumber, const std::string& modelTitle);
    // starts a new chain in the current model, name becomes _atom_site.auth_asym_id
    void addChain(const std::string& chainName);
    // appends a C-alpha only residue to the current chain
    void addCalpha(const char* residueName, int seqId, float x, float y, float z);
    // serializes the accumulated structure as mmCIF into stream, false if gemmi
    // rejected the structure. Stream errors are left to the caller to check.
    bool writeCif(std::ostream& stream);

private:
    GemmiWriter(const GemmiWriter&);
    GemmiWriter& operator=(const GemmiWriter&);
    void* state;
};

#endif //FOLDSEEK_GEMMIWRAPPER_H
