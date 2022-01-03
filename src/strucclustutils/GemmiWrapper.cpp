//
// Created by Martin Steinegger on 6/7/21.
//
#include "GemmiWrapper.h"
#include "mmread.hpp"
#include "gz.hpp"
#include "input.hpp"

//#include "FileUtil.h"

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

bool GemmiWrapper::load(std::string & filename){
    try {
        gemmi::Structure st = openStructure(filename);
        title.clear();
        chain.clear();
        names.clear();
        chainNames.clear();
        ca.clear();
        c.clear();
        cb.clear();
        n.clear();
        ami.clear();
        title.append(  st.get_info("_struct.title"));
        size_t currPos = 0;
        for (gemmi::Model& model : st.models){
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

                    for (gemmi::Atom &atom : res.atoms) {
                        if (atom.name == "CA") {
                            ca_atom.x = atom.pos.x;
                            ca_atom.y = atom.pos.y;
                            ca_atom.z = atom.pos.z;
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
    } catch (std::runtime_error& e) {
        return false;
    }
    return true;
}
