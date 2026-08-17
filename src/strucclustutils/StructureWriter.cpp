#include "StructureWriter.h"

#define GEMMI_WRITE_IMPLEMENTATION
#include "to_mmcif.hpp"
#include "to_cif.hpp"
#include "polyheur.hpp"

#include <cmath>
#include <fstream>
#include <utility>
#include <vector>

namespace {

struct WriterState {
    gemmi::Structure st;
    // model number and source entry title, mmCIF has no per model title item
    std::vector<std::pair<std::string, std::string>> modelTitles;
};

WriterState& state(void* structure) {
    return *static_cast<WriterState*>(structure);
}

// Only the categories that carry information for a C-alpha only structure, the
// default groups would emit empty _exptl, _refine, _cell, ... categories.
// gemmi prints coordinates with the shortest round trip representation of the
// double, so a float 3.8 would end up as 3.79999995. The source files and the
// PDB output carry three decimals, keep the same precision here.
double roundToPdbPrecision(float value) {
    return std::round((double)value * 1000.0) / 1000.0;
}

gemmi::MmcifOutputGroups calphaOutputGroups() {
    gemmi::MmcifOutputGroups groups(false);
    groups.block_name = true;
    groups.entry = true;
    groups.title_keywords = true;
    groups.entity = true;
    groups.struct_asym = true;
    groups.atom_type = true;
    groups.atoms = true;
    groups.group_pdb = true;
    return groups;
}

}

StructureWriter::StructureWriter() : structure(new WriterState()) {}

StructureWriter::~StructureWriter() {
    delete static_cast<WriterState*>(structure);
}

void StructureWriter::reset(const std::string& entryName, const std::string& title) {
    WriterState& s = state(structure);
    s.st = gemmi::Structure();
    s.modelTitles.clear();
    s.st.name = entryName;
    s.st.info["_entry.id"] = entryName;
    if (title.empty() == false) {
        s.st.info["_struct.title"] = title;
    }
}

void StructureWriter::addModel(int modelNumber, const std::string& modelTitle) {
    WriterState& s = state(structure);
    // gemmi keeps the model number as a string
    std::string modelName = std::to_string(modelNumber);
    s.st.models.emplace_back(modelName);
    if (modelTitle.empty() == false) {
        s.modelTitles.emplace_back(modelName, modelTitle);
    }
}

void StructureWriter::addChain(const std::string& chainName) {
    state(structure).st.models.back().chains.emplace_back(chainName);
}

void StructureWriter::addCalpha(const char* residueName, int seqId, float x, float y, float z) {
    gemmi::Residue residue;
    residue.name = residueName;
    residue.seqid = gemmi::SeqId(seqId, ' ');
    residue.label_seq = seqId;
    residue.het_flag = 'A';
    residue.entity_type = gemmi::EntityType::Polymer;

    gemmi::Atom atom;
    atom.name = "CA";
    atom.element = gemmi::El::C;
    atom.occ = 1.0f;
    atom.b_iso = 0.0f;
    atom.pos = gemmi::Position(roundToPdbPrecision(x), roundToPdbPrecision(y), roundToPdbPrecision(z));
    residue.atoms.push_back(atom);

    state(structure).st.models.back().chains.back().residues.push_back(residue);
}

bool StructureWriter::writeCif(const std::string& path) {
    WriterState& s = state(structure);
    // gemmi throws, the callers are compiled without exceptions
    try {
        gemmi::setup_entities(s.st);
        gemmi::cif::Document doc = gemmi::make_mmcif_document(s.st, calphaOutputGroups());
        if (s.modelTitles.empty() == false && doc.blocks.empty() == false) {
            gemmi::cif::Loop& loop = doc.blocks[0].init_loop("_foldseek_model.",
                                                             {"pdbx_PDB_model_num", "title"});
            for (size_t i = 0; i < s.modelTitles.size(); ++i) {
                loop.add_row({gemmi::cif::quote(s.modelTitles[i].first),
                              gemmi::cif::quote(s.modelTitles[i].second)});
            }
        }
        std::ofstream os(path.c_str());
        if (os.is_open() == false) {
            return false;
        }
        gemmi::cif::write_cif_to_stream(os, doc, gemmi::cif::Style::Pdbx);
        os.close();
        return os.fail() == false;
    } catch (const std::exception&) {
        return false;
    }
}
