/**
 * File: amino_acid.h
 * Project: foldcomp
 * Created: 2021-02-04 13:32:06
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Geometric information of amino acids.
 * ---
 * Last Modified: 2022-09-13 15:14:30
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#pragma once
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

class AminoAcid {
public:
    char abb1;
    std::string abb3;
    std::string fullName;
    std::vector<std::string> atoms;
    std::vector<std::string> backboneAtoms = {"N", "CA", "C"};
    std::vector<std::string> sideChainAtoms;
    std::vector<std::string> altAtoms;

    // sideChain uses atom as key and
    std::map < std::string, std::vector<std::string> > sideChain;

    // Geometry from Peptide builder
    std::map<std::string, float> bondLengths;
    std::map<std::string, float> bondAngles;
    std::map<std::string, float> torsionAngles;

    AminoAcid() = default;

    // constructors & destructor
    AminoAcid(char ab1, std::string ab3, std::string name):
        abb1(ab1), abb3(ab3), fullName(name) {}
    AminoAcid(
        char ab1, std::string ab3, std::string name,
        std::vector<std::string> atms,
        std::map< std::string, std::vector<std::string> > sc,
        std::vector<std::string> alt
     ): abb1(ab1), abb3(ab3), fullName(name), atoms(atms), altAtoms(alt), sideChain(sc) {
        for (std::string atm : atms) {
            if (atm != "N" && atm != "CA" && atm != "C") {
                this->sideChainAtoms.push_back(atm);
            }
        }
    };
    AminoAcid(
        char ab1, std::string ab3, std::string name,
        std::vector<std::string> atms,
        std::map< std::string, std::vector<std::string> > sc
    ): abb1(ab1), abb3(ab3), fullName(name), atoms(atms), sideChain(sc) {
        for (std::string atm : atms) {
            if (atm != "N" && atm != "CA" && atm != "C") {
                this->sideChainAtoms.push_back(atm);
            }
        }
    };

    static std::map<std::string, AminoAcid> AminoAcids() {
        std::map<std::string, AminoAcid> output;
        output.emplace("ALA", AminoAcid('A', "ALA", "Alanine", // name
            {"N", "CA", "C", "O", "CB",}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}}}, // sidechain
            {"N", "CA", "C", "CB", "O"}));
        output["ALA"].bondLengths = {{"CA_CB", 1.52}, {"C_O", 1.23}};
        output["ALA"].bondAngles = {{"CA_C_O", 120.31}, {"C_CA_CB", 110.852}};
        // Arginine (R/ARG)
        output.emplace("ARG", AminoAcid(
            'R', "ARG", "Arginine", // name
            { "N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2" }, // atoms
            {{"O", {"N", "CA", "C"}},{"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}},{"CD", {"CA", "CB", "CG"}},
             {"NE", {"CB", "CG", "CD"}},{"CZ", {"CG", "CD", "NE"}},
             {"NH1", {"CD", "NE", "CZ"}},{"NH2", {"CD", "NE", "CZ"}}}, // sidechain
            { "N", "CA", "C", "CB", "O", "CG", "CD", "NE", "NH1", "NH2", "CZ"}
        ));
        output["ARG"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.53},
            {"CG_CD", 1.52}, {"CD_NE", 1.46}, {"NE_CZ", 1.32},
            {"CZ_NH1", 1.31}, {"CZ_NH2", 1.31}
        };
        output["ARG"].bondAngles = {
            {"CA_C_O", 119.745}, {"C_CA_CB", 110.579}, {"CA_CB_CG", 113.233},
            {"CB_CG_CD", 110.787}, {"CG_CD_NE", 111.919}, {"CD_NE_CZ", 125.192},
            {"NE_CZ_NH1", 120.077}, {"NE_CZ_NH2", 120.077}
        };
        // Asparagine (N/ASN)
        output.emplace("ASN", AminoAcid(
            'N', "ASN", "Asparagine", // name
            { "N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"}, // atoms
            {{"O", {"N", "CA", "C"}},{"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}},{"OD1", {"CA", "CB", "CG"}},
             {"ND2", {"CA", "CB", "CG"}}}, // sidechain
            { "N", "CA", "C", "CB", "O", "CG", "ND2", "OD1"}
        ));
        output["ASN"].bondLengths = {
            {"CA_CB", 1.52}, {"C_O", 1.23}, {"CB_CG", 1.52}, {"CG_OD1", 1.23}, {"CG_ND2", 1.325}
        };
        output["ASN"].bondAngles = {
            {"CA_C_O", 120.313}, {"C_CA_CB", 110.852}, {"CA_CB_CG", 113.232},
            {"CB_CG_OD1", 120.85}, {"CB_CG_ND2", 116.48}
        };
        // Aspartic Acid (D/ASP)
        output.emplace("ASP", AminoAcid(
            'D', "ASP", "Aspartic acid", // name
            { "N", "CA", "C", "O", "CB", "CG", "OD1", "OD2" }, // atoms
            {{"O", {"N", "CA", "C"}},{"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}},{"OD1", {"CA", "CB", "CG"}},
             {"OD2", {"CA", "CB", "CG"}}}, // sidechain
            { "N", "CA", "C", "CB", "O", "CG", "OD1", "OD2"}
        ));
        output["ASP"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.52}, {"CG_OD1", 1.248}, {"CG_OD2", 1.248}
        };
        output["ASP"].bondAngles = {
            {"CA_C_O", 121.051}, {"C_CA_CB", 110.871}, {"CA_CB_CG", 113.232},
            {"CB_CG_OD1", 118.344}, {"CB_CG_OD2", 118.344}
        };
        // Cysteine (C/CYS)
        output.emplace("CYS", AminoAcid(
            'C', "CYS", "Cysteine", // name
            {"N", "CA", "C", "O", "CB", "SG"}, // atoms
            {{"O", {"N", "CA", "C"}},{"CB", {"O", "C", "CA"}},
             {"SG", {"N", "CA", "CB"}}}, // sidechain
            { "N", "CA", "C", "CB", "O", "SG"}
        ));
        output["CYS"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_SG", 1.8}
        };
        output["CYS"].bondAngles = {
            {"CA_C_O", 120.063}, {"C_CA_CB", 111.078}, {"CA_CB_SG", 113.817}
        };
        // Glutamine (Q/GLN)
        output.emplace("GLN", AminoAcid(
            'Q', "GLN", "Glutamine", // name
            {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"}, // atoms
            {{"O", {"N", "CA", "C"}},{"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}},{"CD", {"CA", "CB", "CG"}},
             {"OE1", {"CB", "CG", "CD"}},{"NE2", {"CB", "CG", "CD"}}}, // sidechain
            { "N", "CA", "C", "CB", "O", "CG", "CD", "NE2", "OE1"}
        ));
        output["GLN"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.52}, {"CG_CD", 1.52},
            {"CD_OE1", 1.23}, {"CD_NE2", 1.32}
        };
        output["GLN"].bondAngles = {
            {"CA_C_O", 120.211}, {"C_CA_CB", 109.5}, {"CA_CB_CG", 113.292},
            {"CB_CG_CD", 112.811}, {"CG_CD_OE1", 121.844}, {"CG_CD_NE2", 116.50}
        };
        // Glutamic Acid (E/GLU)
        output.emplace("GLU", AminoAcid(
            'E', "GLU", "Glutamic acid",
            {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
            {{"O", {"N", "CA", "C"}},{"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}},{"CD", {"CA", "CB", "CG"}},
             {"OE1", {"CB", "CG", "CD"}},{"OE2", {"CB", "CG", "CD"}}},
            { "N", "CA", "C", "CB", "O", "CG", "CD", "OE1", "OE2"}
        ));
        output["GLU"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.52}, {"CG_CD", 1.52},
            {"CD_OE1", 1.25}, {"CD_OE2", 1.25}
        };
        output["GLU"].bondAngles = {
            {"CA_C_O", 120.594}, {"C_CA_CB", 110.538}, {"CA_CB_CG", 113.82},
            {"CB_CG_CD", 112.912}, {"CG_CD_OE1", 118.479}, {"CG_CD_OE2", 118.479}
        };
        // Glycine (G/GLY)
        output.emplace("GLY", AminoAcid(
            'G', "GLY", "Glycine", // name
            {"N", "CA", "C", "O"}, // atoms
            {{"O", {"N", "CA", "C"}}}, // sidechain
            { "N", "CA", "C", "O"}
        ));
        output["GLY"].bondLengths = {{"C_O", 1.23}};
        output["GLY"].bondAngles = {{"CA_C_O", 120.522}};
        // Histidine (H/HIS)
        output.emplace("HIS", AminoAcid(
            'H', "HIS", "Histidine",
            { "N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2" }, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"ND1", {"CA", "CB", "CG"}},
             {"CD2", {"CA", "CB", "CG"}}, {"CE1", {"CB", "CG", "ND1"}},
             {"NE2", {"CB", "CG", "CD2"}}},
            { "N", "CA", "C", "CB", "O", "CG", "CD2", "ND1", "CE1", "NE2" }
        ));
        output["HIS"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.5}, {"CG_ND1", 1.38},
            {"CG_CD2", 1.36}, {"ND1_CE1", 1.33}, {"CD2_NE2", 1.38}
        };
        output["HIS"].bondAngles = {
            {"CA_C_O", 120.548}, {"C_CA_CB", 111.329}, {"CA_CB_CG", 113.468},
            {"CB_CG_CD2", 130.61}, {"CB_CG_ND1", 122.85}, {"CG_CD2_NE2", 107.439},
            {"CG_ND1_CE1", 108.589}
        };
        // Isoleucine (I/ILE)
        output.emplace("ILE", AminoAcid(
            'I', "ILE", "Isoleucine", // name
            {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG1", {"N", "CA", "CB"}}, {"CG2", {"N", "CA", "CB"}},
             {"CD1", {"CA", "CB", "CG1"}}},
             { "N", "CA", "C", "CB", "O", "CG1", "CG2", "CD1" }
        ));
        output["ILE"].bondLengths = {
            {"CA_CB", 1.54}, {"C_O", 1.235}, {"CB_CG1", 1.53}, {"CB_CG2", 1.52},
            {"CG1_CD1", 1.51}
        };
        output["ILE"].bondAngles = {
            {"CA_C_O", 120.393}, {"C_CA_CB", 111.983}, {"CA_CB_CG1", 110.5},
            {"CA_CB_CG2", 110.5}, {"CB_CG1_CD1", 113.97}
        };
        // Leucine (L/LEU)
        output.emplace("LEU", AminoAcid(
            'L', "LEU", "Leucine",
            { "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2" }, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"CD1", {"CA", "CB", "CG"}},
             {"CD2", {"CA", "CB", "CG"}} },
            { "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2" }
        ));
        output["LEU"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.235}, {"CB_CG", 1.53}, {"CG_CD1", 1.52},
            {"CG_CD2", 1.52}
        };
        output["LEU"].bondAngles = {
            {"CA_C_O", 120.211}, {"C_CA_CB", 110.418}, {"CA_CB_CG", 116.10},
            {"CB_CG_CD1", 110.58}, {"CB_CG_CD2", 110.58}
        };
        // Lysine (K/LYS)
        // 2022-06-10 21:56:59 - TODO: RECALCULATE GEOMETRY FOR LYSINE
        output.emplace("LYS", AminoAcid(
            'K', "LYS", "Lysine",
            { "N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ" }, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"CD", {"CA", "CB", "CG"}},
             {"CE", {"CB", "CG", "CD"}}, {"NZ", {"CG", "CD", "CE"}}},
            { "N", "CA", "C", "CB", "O", "CG", "CD", "CE", "NZ" }
        ));
        output["LYS"].bondLengths = {
            {"C_O", 1.23}, {"CA_CB", 1.53}, {"CB_CG", 1.52}, {"CG_CD", 1.52},
            {"CD_CE", 1.52}, {"CE_NZ", 1.49} // sidechain
        };
        output["LYS"].bondAngles = {
            {"CA_C_O", 120.54}, {"C_CA_CB", 109.5}, {"CA_CB_CG", 113.83},
            {"CB_CG_CD", 111.79}, {"CG_CD_CE", 111.79}, {"CD_CE_NZ", 112.25}
        };
        // Methionine (M/MET)
        output.emplace("MET", AminoAcid(
            'M', "MET", "Methionine",
            { "N", "CA", "C", "O", "CB", "CG", "SD", "CE"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"SD", {"CA", "CB", "CG"}},
             {"CE", {"CB", "CG", "SD"}}},
            { "N", "CA", "C", "CB", "O", "CG", "SD", "CE" }
        ));
        output["MET"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.52}, {"CG_SD", 1.8},
            {"SD_CE", 1.79}
        };
        output["MET"].bondAngles = {
            {"CA_C_O", 120.148}, {"C_CA_CB", 110.833}, {"CA_CB_CG", 113.68},
            {"CB_CG_SD", 112.773}, {"CG_SD_CE", 100.61}
        };
        // Phenylalanine (F/PHE)
        output.emplace("PHE", AminoAcid(
            'F', "PHE", "Phenylalanine",
            { "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"CD1", {"CA", "CB", "CG"}},
             {"CD2", {"CA", "CB", "CG"}}, {"CE1", {"CB", "CG", "CD1"}},
             {"CE2", {"CB", "CG", "CD2"}}, {"CZ", {"CG", "CD1", "CE1"}}},
            { "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2", "CE1", "CE2", "CZ" }
        ));
        output["PHE"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.51}, {"CG_CD1", 1.385},
            {"CG_CD2", 1.385}, {"CD1_CE1", 1.385}, {"CD2_CE2", 1.385},
            {"CE1_CZ", 1.385}
        };
        output["PHE"].bondAngles = {
            {"CA_C_O", 120.283}, {"C_CA_CB", 110.846}, {"CA_CB_CG", 114.0},
            {"CB_CG_CD1", 120.0}, {"CB_CG_CD2", 120.0}, {"CG_CD1_CE1", 120.0},
            {"CG_CD2_CE2", 120.0}, {"CD1_CE1_CZ", 120.0}
        };
        // Proline (P/PRO)
        output.emplace("PRO", AminoAcid(
            'P', "PRO", "Proline", // name
            { "N", "CA", "C", "O", "CB", "CG", "CD"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"CD", {"CA", "CB", "CG"}}},
            { "N", "CA", "C", "CB", "O", "CG", "CD" } // atoms
        ));
        output["PRO"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.49}, {"CG_CD", 1.50}
        };
        output["PRO"].bondAngles = {
            {"CA_C_O", 120.6}, {"C_CA_CB", 111.372}, {"CA_CB_CG", 104.21},
            {"CB_CG_CD", 105.0}
        };
        // Serine (S/SER)
        output.emplace("SER", AminoAcid(
            'S', "SER", "Serine",
            { "N", "CA", "C", "O", "CB", "OG"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"OG", {"N", "CA", "CB"}}},
            { "N", "CA", "C", "CB", "O", "OG" }
        ));
        output["SER"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_OG", 1.417}
        };
        output["SER"].bondAngles = {
            {"CA_C_O", 120.475}, {"C_CA_CB", 110.248}, {"CA_CB_OG", 111.132}
        };
        // Threonine (T/THR)
        output.emplace("THR", AminoAcid(
            'T', "THR", "Threonine",
            {"N", "CA", "C", "O", "CB", "OG1", "CG2"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"OG1", {"N", "CA", "CB"}}, {"CG2", {"N", "CA", "CB"}}},
            { "N", "CA", "C", "CB", "O", "CG2", "OG1" }
        ));
        output["THR"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_OG1", 1.43}, {"CB_CG2", 1.52}
        };
        output["THR"].bondAngles = {
            {"CA_C_O", 120.252}, {"C_CA_CB", 110.075}, {"CA_CB_OG1", 109.442},
            {"CA_CB_CG2", 111.457}
        };
        // Tryptophan (W/TRP)
        output.emplace("TRP", AminoAcid(
            'W', "TRP", "Tryptophan",
            {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2",
             "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"CD1", {"CA", "CB", "CG"}},
             {"CD2", {"CA", "CB", "CG"}}, {"NE1", {"CB", "CG", "CD1"}},
             {"CE2", {"CB", "CG", "CD2"}}, {"CE3", {"CB", "CG", "CD2"}},
             {"CZ2", {"CG", "CD2", "CE2"}}, {"CZ3", {"CG", "CD2", "CE3"}},
             {"CH2", {"CD2", "CE2", "CZ2"}}},
            { "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2",
              "CE2", "CE3", "NE1", "CH2", "CZ2", "CZ3" }
        ));
        output["TRP"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.23}, {"CB_CG", 1.50},
            {"CG_CD1", 1.36}, {"CG_CD2", 1.44}, {"CD1_NE1", 1.38},
            {"CD2_CE2", 1.41}, {"CD2_CE3", 1.40}, {"CE2_CZ2", 1.40},
            {"CE3_CZ3", 1.384}, {"CZ2_CH2", 1.367}
        };
        output["TRP"].bondAngles = {
            {"CA_C_O", 120.178}, {"C_CA_CB", 110.852}, {"CA_CB_CG", 114.10},
            {"CB_CG_CD1", 126.712}, {"CB_CG_CD2", 126.712}, {"CG_CD1_NE1", 109.959},
            {"CG_CD2_CE2", 107.842}, {"CG_CD2_CE3", 133.975}, {"CD2_CE2_CZ2", 120.0},
            {"CD2_CE3_CZ3", 120.0}, {"CE2_CZ2_CH2", 120.0}
        };
        // Tyrosine (Y/TYR)
        output.emplace("TYR", AminoAcid(
            'Y', "TYR", "Tyrosine", // name
            {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2",
             "CE1", "CE2", "CZ", "OH"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG", {"N", "CA", "CB"}}, {"CD1", {"CA", "CB", "CG"}},
             {"CD2", {"CA", "CB", "CG"}}, {"CE1", {"CB", "CG", "CD1"}},
             {"CE2", {"CB", "CG", "CD2"}}, {"CZ", {"CG", "CD1", "CE1"}},
             {"OH", {"CD1", "CE1", "CZ"}}},
            { "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2",
              "CE1", "CE2", "OH", "CZ" }
        ));
        output["TYR"].bondLengths = {
            {"CA_CB", 1.53}, {"C_O", 1.235}, {"CB_CG", 1.51},
            {"CG_CD1", 1.39}, {"CG_CD2", 1.39}, {"CD1_CE1", 1.38},
            {"CD2_CE2", 1.38}, {"CE1_CZ", 1.378}, {"CZ_OH", 1.375}
        };
        output["TYR"].bondAngles = {
            {"CA_C_O", 120.608}, {"C_CA_CB", 110.852}, {"CA_CB_CG", 113.744},
            {"CB_CG_CD1", 120.937}, {"CB_CG_CD2", 120.937}, {"CG_CD1_CE1", 120.0},
            {"CG_CD2_CE2", 120.0}, {"CD1_CE1_CZ", 120.0}, {"CE1_CZ_OH", 120.0}
        };
        // Valine (V/VAL)
        output.emplace("VAL", AminoAcid(
            'V', "VAL", "Valine", // name
            { "N", "CA", "C", "O", "CB", "CG1", "CG2"}, // atoms
            {{"O", {"N", "CA", "C"}}, {"CB", {"O", "C", "CA"}},
             {"CG1", {"N", "CA", "CB"}}, {"CG2", {"N", "CA", "CB"}}},
            { "N", "CA", "C", "CB", "O", "CG1", "CG2" }
        ));
        output["VAL"].bondLengths = {
            {"CA_CB", 1.54}, {"C_O", 1.235}, {"CB_CG1", 1.52}, {"CB_CG2", 1.52}
        };
        output["VAL"].bondAngles = {
            {"CA_C_O", 120.472}, {"C_CA_CB", 111.381},
            {"CA_CB_CG1", 110.7}, {"CA_CB_CG2", 110.4}
        };
        // output.emplace("ASX", AminoAcid('B', "ASX", "Asparagine/aspartic acid"));
        // output.emplace("GLX", AminoAcid('Z', "GLX", "Glutamine/glutamic acid"));
        output.emplace("UNK", AminoAcid('X', "UNK", "Unknown"));
        return output;
    }

    std::map<unsigned int, std::string> AminoAcidIndexMap() {
        std::map<unsigned int, std::string> output;
        output[0]="ALA";
        output[1]="ARG";
        output[2]="ASN";
        output[3]="ASP";
        output[4]="CYS";
        output[5]="GLN";
        output[6]="GLU";
        output[7]="GLY";
        output[8]="HIS";
        output[9]="ILE";
        output[10]="LEU";
        output[11]="LYS";
        output[12]="MET";
        output[13]="PHE";
        output[14]="PRO";
        output[15]="SER";
        output[16]="THR";
        output[17]="TRP";
        output[18]="TYR";
        output[19]="VAL";
        // output[20]="ASX";
        // output[21]="GLX";
        return output;
    }
    std::map<std::string, unsigned int> IndexAminoAcidMap() {
        std::map<std::string, unsigned int> output;
        std::map<unsigned int, std::string> aa_ind_map = this->AminoAcidIndexMap();
        std::map<unsigned int, std::string>::iterator it;
        for (it = aa_ind_map.begin(); it != aa_ind_map.end(); it++) {
            output[it->second] = it->first;
        }
        return output;
    }

};

int writeAminoAcidMapToFile(std::map<std::string, AminoAcid>& aa_map, std::string filename);

std::vector<std::string> getAminoAcidList(void);