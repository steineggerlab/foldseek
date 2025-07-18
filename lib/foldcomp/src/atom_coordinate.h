/**
 * File: atom_coordinate.h
 * Project: foldcomp
 * Created: 2021-01-18 12:43:08
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     The data type to handle atom coordinate comes here.
 * ---
 * Last Modified: 2022-11-29 14:39:06
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include "float3d.h"
#include "tcbspan.h"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

class AtomCoordinate {
public:
    AtomCoordinate() = default;
    AtomCoordinate(
        std::string a, std::string r, std::string c,
        int ai, int ri, float x, float y, float z,
        float occupancy = 0.0f, float tempFactor = 0.0f
    );
    AtomCoordinate(
        std::string a, std::string r, std::string c,
        int ai, int ri, float3d coord,
        float occupancy = 0.0f, float tempFactor = 0.0f
    );
    // data
    std::string atom;
    std::string residue;
    std::string chain;
    int atom_index;
    int residue_index;
    float3d coordinate;
    float occupancy;
    float tempFactor;

    // operators
    bool operator==(const AtomCoordinate& other) const;
    bool operator!=(const AtomCoordinate& other) const;
    //BackboneChain toCompressedResidue();

    //methods
    bool isBackbone() const;
    void print(int option = 0) const ;
    void setTempFactor(float tf) { this->tempFactor = tf; };
};

std::vector<float3d> extractCoordinates(const std::vector<AtomCoordinate>& atoms);

static inline void extractCoordinates(
    float3d* output,
    const AtomCoordinate& atom1,
    const AtomCoordinate& atom2,
    const AtomCoordinate& atom3
) {
    output[0] = atom1.coordinate;
    output[1] = atom2.coordinate;
    output[2] = atom3.coordinate;
}

std::vector<AtomCoordinate> extractChain(
    std::vector<AtomCoordinate>& atoms, std::string chain
);

std::vector<AtomCoordinate> filterBackbone(const tcb::span<AtomCoordinate>& atoms);

void printAtomCoordinateVector(std::vector<AtomCoordinate>& atoms, int option = 0);

std::vector<AtomCoordinate> weightedAverage(
    const std::vector<AtomCoordinate>& origAtoms, const std::vector<AtomCoordinate>& revAtoms
);

void writeAtomCoordinatesToPDB(
    std::vector<AtomCoordinate>& atoms, std::string title, std::ostream& pdb_path
);
int writeAtomCoordinatesToPDBFile(
    std::vector<AtomCoordinate>& atoms, std::string title, std::string pdb_path
);

std::vector<std::vector<AtomCoordinate>> splitAtomByResidue(
    const tcb::span<AtomCoordinate>& atomCoordinates
);

std::vector<std::string> getResidueNameVector(
    const tcb::span<AtomCoordinate>& atomCoordinates
);

AtomCoordinate findFirstAtom(const std::vector<AtomCoordinate>& atoms, std::string atom_name);
void setAtomIndexSequentially(std::vector<AtomCoordinate>& atoms, int start);
void removeAlternativePosition(std::vector<AtomCoordinate>& atoms);

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, int residue_index,
    std::vector<std::string> atomNames = {"N", "CA", "C"}
);

std::vector<std::vector<AtomCoordinate>> getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, std::vector<int> residue_index,
    std::vector<std::string> atomNames = {"N", "CA", "C"}
);
float RMSD(std::vector<AtomCoordinate>& atoms1, std::vector<AtomCoordinate>& atoms2);

template <int32_t T, int32_t P>
void ftoa(float n, char* s);

std::vector<std::pair<size_t, size_t>> identifyChains(const std::vector<AtomCoordinate>& atoms);
std::vector<std::pair<size_t, size_t>> identifyDiscontinousResInd(const std::vector<AtomCoordinate>& atoms, size_t chain_start, size_t chain_end);
