//
// Created by Martin Steinegger on 6/7/21.
//
#include "GemmiWrapper.h"
#include "mmread.hpp"
#include "gz.hpp"
#include "input.hpp"
#include "foldcomp.h"

GemmiWrapper::GemmiWrapper(){
    threeAA2oneAA = {{"ALA",'A'},  {"ARG",'R'},  {"ASN",'N'}, {"ASP",'D'},
                     {"CYS",'C'},  {"GLN",'Q'},  {"GLU",'E'}, {"GLY",'G'},
                     {"HIS",'H'},  {"ILE",'I'},  {"LEU",'L'}, {"LYS",'K'},
                     {"MET",'M'},  {"PHE",'F'},  {"PRO",'P'}, {"SER",'S'},
                     {"THR",'T'},  {"TRP",'W'},  {"TYR",'Y'}, {"VAL",'V'},
                     {"UNK",'X'}};
}

gemmi::Structure openStructure(const std::string & filename){
    gemmi::MaybeGzipped infile(filename);
    gemmi::CoorFormat format = gemmi::coor_format_from_ext(infile.basepath());
    if(format != gemmi::CoorFormat::Unknown && format != gemmi::CoorFormat::Unknown){
        return gemmi::read_structure(infile, format);
    }else{
        return gemmi::read_structure(infile, gemmi::CoorFormat::Pdb);
    }
}

// https://stackoverflow.com/questions/1448467/initializing-a-c-stdistringstream-from-an-in-memory-buffer/1449527
struct OneShotReadBuf : public std::streambuf
{
    OneShotReadBuf(char* s, std::size_t n)
    {
        setg(s, s, s + n);
    }
};

bool GemmiWrapper::loadFromBuffer(const char * buffer, size_t bufferSize, std::string & name) {
    if (gemmi::iends_with(name, ".fcz")) {
        OneShotReadBuf buf((char *) buffer, bufferSize);
        std::istream istr(&buf);
        if (!istr) {
            return false;
        }
        return loadFoldcompStructure(istr);
    }
    try {
        gemmi::MaybeGzipped infile(name);
        gemmi::CoorFormat format = gemmi::coor_format_from_ext(infile.basepath());
        gemmi::Structure st;
        switch (format) {
            case gemmi::CoorFormat::Pdb:
                st = gemmi::pdb_impl::read_pdb_from_stream(gemmi::MemoryStream(buffer, bufferSize), name, gemmi::PdbReadOptions());
                break;
            case gemmi::CoorFormat::Mmcif:
                st = gemmi::make_structure(gemmi::cif::read_memory(buffer, bufferSize, name.c_str()));
                break;
            case gemmi::CoorFormat::Unknown:
            case gemmi::CoorFormat::Detect:
                return false;
        }
        updateStructure((void*) &st, name);
    } catch (std::runtime_error& e) {
        return false;
    }
    return true;
}

bool GemmiWrapper::loadFoldcompStructure(std::istream& stream) {
    std::cout.setstate(std::ios_base::failbit);
    Foldcomp fc;
    int res = fc.read(stream);
    if (res != 0) {
        return false;
    }
    std::vector<AtomCoordinate> coordinates;
    fc.useAltAtomOrder = false;
    res = fc.decompress(coordinates);
    if (res != 0) {
        return false;
    }
    std::cout.clear();
    if (coordinates.size() == 0) {
        return false;
    }

    title.clear();
    chain.clear();
    names.clear();
    chainNames.clear();
    ca.clear();
    ca_bfactor.clear();
    c.clear();
    cb.clear();
    n.clear();
    ami.clear();
    title.append(fc.strTitle);
    names.push_back(fc.strTitle);
    const AtomCoordinate& first = coordinates[0];
    chainNames.push_back(first.chain);
    int residueIndex = INT_MAX;
    Vec3 ca_atom = {NAN, NAN, NAN};
    Vec3 cb_atom = {NAN, NAN, NAN};
    Vec3 n_atom  = {NAN, NAN, NAN};
    Vec3 c_atom  = {NAN, NAN, NAN};
    float ca_atom_bfactor = 0.0;
    for (std::vector<AtomCoordinate>::const_iterator it = coordinates.begin(); it != coordinates.end(); ++it) {
        const AtomCoordinate& atom = *it;
        if (atom.residue_index != residueIndex) {
            if (residueIndex != INT_MAX) {
                ca.push_back(ca_atom);
                cb.push_back(cb_atom);
                n.push_back(n_atom);
                c.push_back(c_atom);
                ca_bfactor.push_back(ca_atom_bfactor);
                ca_atom = {NAN, NAN, NAN};
                cb_atom = {NAN, NAN, NAN};
                n_atom  = {NAN, NAN, NAN};
                c_atom  = {NAN, NAN, NAN};
                ca_atom_bfactor = 0.0;
            }
            if (threeAA2oneAA.find(atom.residue) == threeAA2oneAA.end()) {
                ami.push_back('X');
            } else {
                ami.push_back(threeAA2oneAA[atom.residue]);
            }
            residueIndex = atom.residue_index;
        }

        if (atom.atom == "CA") {
            ca_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
            ca_atom_bfactor = atom.tempFactor;
        } else if (atom.atom == "CB") {
            cb_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
        } else if (atom.atom == "N") {
            n_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
        } else if (atom.atom == "C") {
            c_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
        }
    }
    ca.push_back(ca_atom);
    cb.push_back(cb_atom);
    n.push_back(n_atom);
    c.push_back(c_atom);
    ca_bfactor.push_back(ca_atom_bfactor);
    chain.emplace_back(0, ca.size());
    return true;
}

void GemmiWrapper::updateStructure(void * void_st, std::string & filename) {
    gemmi::Structure * st = (gemmi::Structure *) void_st;

    title.clear();
    chain.clear();
    names.clear();
    chainNames.clear();
    ca.clear();
    ca_bfactor.clear();
    c.clear();
    cb.clear();
    n.clear();
    ami.clear();
    title.append(  st->get_info("_struct.title"));
    size_t currPos = 0;
    for (gemmi::Model& model : st->models){
        for (gemmi::Chain& ch : model.chains) {

            size_t chainStartPos = currPos;
            size_t pos = filename.find_last_of("\\/");
            std::string name = (std::string::npos == pos)
                               ? filename
                               : filename.substr(pos+1, filename.length());
            //name.push_back('_');
            chainNames.push_back(ch.name);
            names.push_back(name);
            for (gemmi::Residue &res : ch.residues) {
                if (res.het_flag != 'A')
                    continue;
                Vec3 ca_atom = {NAN, NAN, NAN};
                Vec3 cb_atom = {NAN, NAN, NAN};
                Vec3 n_atom  = {NAN, NAN, NAN};
                Vec3 c_atom  = {NAN, NAN, NAN};
                float ca_atom_bfactor;
                for (gemmi::Atom &atom : res.atoms) {
                    if (atom.name == "CA") {
                        ca_atom.x = atom.pos.x;
                        ca_atom.y = atom.pos.y;
                        ca_atom.z = atom.pos.z;
                        ca_atom_bfactor = atom.b_iso;
                    } else if (atom.name == "CB") {
                        cb_atom.x = atom.pos.x;
                        cb_atom.y = atom.pos.y;
                        cb_atom.z = atom.pos.z;
                    } else if (atom.name == "N") {
                        n_atom.x = atom.pos.x;
                        n_atom.y = atom.pos.y;
                        n_atom.z = atom.pos.z;
                    } else if (atom.name == "C") {
                        c_atom.x = atom.pos.x;
                        c_atom.y = atom.pos.y;
                        c_atom.z = atom.pos.z;
                    }
                }
                ca_bfactor.push_back(ca_atom_bfactor);
                ca.push_back(ca_atom);
                cb.push_back(cb_atom);
                n.push_back(n_atom);
                c.push_back(c_atom);
                currPos++;
                if (threeAA2oneAA.find(res.name) == threeAA2oneAA.end()) {
                    ami.push_back('X');
                } else {
                    ami.push_back(threeAA2oneAA[res.name]);
                }
            }
            chain.push_back(std::make_pair(chainStartPos, currPos));
        }
    }
}

bool GemmiWrapper::load(std::string & filename){
    if (gemmi::iends_with(filename, ".fcz")) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            return false;
        }
        return loadFoldcompStructure(in);
    }
    try {
        gemmi::Structure st = openStructure(filename);
        updateStructure((void*) &st, filename);
    } catch (std::runtime_error& e) {
        return false;
    }
    return true;
}
