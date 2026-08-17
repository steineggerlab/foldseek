#ifndef STRUCTUREWRITER_H
#define STRUCTUREWRITER_H

#include <cstddef>
#include <string>

// Writes structures with gemmi. The gemmi headers need exceptions, which the
// foldseek framework is compiled without, so the implementation lives in the
// gemmiwrapper library and only this plain interface is exposed here.
//
// Usage: reset(), then addModel()/addChain()/addCalpha() in that order,
// then writeCif(). A writer can be reused after writeCif().
class StructureWriter {
public:
    StructureWriter();
    ~StructureWriter();

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
    // serializes the accumulated structure as mmCIF, false if it could not be written
    bool writeCif(const std::string& path);

private:
    StructureWriter(const StructureWriter&);
    StructureWriter& operator=(const StructureWriter&);
    void* structure;
};

#endif
