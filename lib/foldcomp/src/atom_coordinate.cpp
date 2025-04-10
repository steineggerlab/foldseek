/**
 * File: atom_coordinate.cpp
 * Project: foldcomp
 * Created: 2021-01-18 12:53:34
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     The data type to handle atom coordinate comes here.
 * ---
 * Last Modified: Fri Mar 03 2023
 * Modified By: Hyunbin Kim
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "atom_coordinate.h"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream> // IWYU pragma: keep

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param x A float for x coordinate
 * @param y A float for y coordinate
 * @param z A float for z coordinate
 */
AtomCoordinate::AtomCoordinate(
    std::string a, std::string r, std::string c,
    int ai, int ri, float x, float y, float z,
    float occupancy, float tempFactor
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri), occupancy(occupancy), tempFactor(tempFactor) {
    this->coordinate = {x, y, z};
}

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param coord A float vector for x,y,z coordinates.
 */
AtomCoordinate::AtomCoordinate(
    std::string a, std::string r, std::string c,
    int ai, int ri, float3d coord,
    float occupancy, float tempFactor
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri), coordinate(coord), occupancy(occupancy), tempFactor(tempFactor) {
}

bool AtomCoordinate::operator==(const AtomCoordinate& other) const {
    return (
        (this->atom == other.atom) &&
        (this->atom_index == other.atom_index) &&
        (this->residue == other.residue) &&
        (this->residue_index == other.residue_index) &&
        (this->chain == other.chain) &&
        (this->coordinate.x == other.coordinate.x) &&
        (this->coordinate.y == other.coordinate.y) &&
        (this->coordinate.z == other.coordinate.z)
    );
}
bool AtomCoordinate::operator!=(const AtomCoordinate& other) const {
    return !(*this == other);
}

bool AtomCoordinate::isBackbone() const {
    return ((this->atom == "N") ||(this->atom == "CA") ||(this->atom == "C"));
}

void AtomCoordinate::print(int option) const {
    std::cout << "Atom: " << this->atom << std::endl;
    if (option != 0) {
        std::cout << "Residue: " << this->residue << std::endl;
        std::cout << "Chain: " << this->chain << std::endl;
        std::cout << "Atom Index: " << this->atom_index << std::endl;
        std::cout << "Residue Index: " << this->residue_index << std::endl;
        if (option == 2) {
            std::cout << "Coordinate: ";
            std::cout << this->coordinate.x << " ";
            std::cout << this->coordinate.y << " ";
            std::cout << this->coordinate.z << " ";
            std::cout << std::endl;
        }
    }
}

/**
 * @brief Extracts coordinates from AtomCoordinate vector
 *
 * @param atoms A vector of AtomCoordinate
 * @return std::vector< std::vector<float> >
 */
std::vector<float3d> extractCoordinates(
    const std::vector<AtomCoordinate>& atoms
) {
    std::vector<float3d> output(atoms.size());
    for (size_t i = 0; i < atoms.size(); i++) {
        output[i] = atoms[i].coordinate;
    }
    return output;
}

std::vector<AtomCoordinate> extractChain(
    std::vector<AtomCoordinate>& atoms, std::string chain
) {
    std::vector<AtomCoordinate> output;
    int total = atoms.size();
    for (int i = 0; i < total; i++) {
        const AtomCoordinate& curr_atm = atoms[i];
        if (i < (total-1)) {
            const AtomCoordinate& next_atm = atoms[i + 1];
            if (next_atm.atom == curr_atm.atom) {
                continue;
            }
        }
        if (curr_atm.chain == chain) {
            output.push_back(curr_atm);
        }
    }
    return output;
}

void printAtomCoordinateVector(std::vector<AtomCoordinate>& atoms, int option) {
    for (const AtomCoordinate& curr_atm : atoms) {
        curr_atm.print(option);
    }
}

std::vector<AtomCoordinate> filterBackbone(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.isBackbone()) {
            output.emplace_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> weightedAverage(
    const std::vector<AtomCoordinate>& origAtoms, const std::vector<AtomCoordinate>& revAtoms
) {
    std::vector<AtomCoordinate> output;
    output.reserve(origAtoms.size());
    int total = origAtoms.size();
    for (int i = 0; i < total; i++) {
        const AtomCoordinate& curr_atm = origAtoms[i];
        const AtomCoordinate&  rev_atm = revAtoms[i];
        output.emplace_back(
            curr_atm.atom, curr_atm.residue, curr_atm.chain,
            curr_atm.atom_index, curr_atm.residue_index,
            ((curr_atm.coordinate.x * (float)(total - i)) + (rev_atm.coordinate.x * (float)i)) / (float)total,
            ((curr_atm.coordinate.y * (float)(total - i)) + (rev_atm.coordinate.y * (float)i)) / (float)total,
            ((curr_atm.coordinate.z * (float)(total - i)) + (rev_atm.coordinate.z * (float)i)) / (float)total
        );
    }
    return output;
}

void reverse(char* s) {
    for (int i = 0, j = strlen(s)-1; i < j; i++, j--) {
        char c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

void itoa_pos_only(int n, char* s) {
    int i = 0;
    do {
        // generate digits in reverse order
        // get next digit
        s[i++] = n % 10 + '0';
    // shift to next
    } while ((n /= 10) > 0);
    s[i] = '\0';
    reverse(s);
}

template <int32_t T, int32_t P>
void fast_ftoa(float n, char* s) {
    float rounded = n + ((n < 0) ? -(0.5f / T) : (0.5f / T));
    int32_t integer = (int32_t)rounded;
    int32_t decimal = (int32_t)((rounded - (float)integer) * (float)T);
    char* data = s;
    if (n < 0) {
        integer = std::abs(integer);
        decimal = std::abs(decimal);
        *data = '-';
        data++;
    }
    itoa_pos_only(integer, data);
    data += strlen(data);
    *data = '.';
    data++;
    char buffer[10];
    itoa_pos_only(decimal, buffer);
    int32_t len = strlen(buffer);
    for (int32_t i = 0; i < (P - len); i++) {
        *data = '0';
        data++;
    }
    memcpy(data, buffer, len);
    // add a null terminator
    data += len;
    *data = '\0';
    // std::string check(s);
    // std::ostringstream ss;
    // ss << std::fixed << std::setprecision(P) << n;
    // if (ss.str() != check) {
    //     std::cout << "ERROR: " << ss.str() << " != " << check << " ORIG: " << std::fixed << std::setprecision(10) << n << std::endl;
    // }
}

void writeAtomCoordinatesToPDB(
    std::vector<AtomCoordinate>& atoms, std::string title, std::ostream& pdb_stream
) {
    // Write title
    // Check if title is too long and if so, write the title in multiple lines
    if (title != "") {
        const char* headerData = title.c_str();
        size_t headerLen = title.length();
        int remainingHeader = headerLen;
        char buffer[128];
        int written = snprintf(buffer, sizeof(buffer), "TITLE     %.*s\n",  std::min(70, (int)remainingHeader), headerData);
        if (written >= 0 && written < (int)sizeof(buffer)) {
            pdb_stream << buffer;
        }
        remainingHeader -= 70;
        int continuation = 2;
        while (remainingHeader > 0) {
            written = snprintf(buffer, sizeof(buffer), "TITLE  % 3d%.*s\n", continuation, std::min(70, (int)remainingHeader), headerData + (headerLen - remainingHeader));
            if (written >= 0 && written < (int)sizeof(buffer)) {
                pdb_stream << buffer;
            }
            remainingHeader -= 70;
            continuation++;
        }
    }

    int total = atoms.size();
    std::string residue;
    for (int i = 0; i < total; i++) {
        pdb_stream << "ATOM  "; // 1-4 ATOM
        pdb_stream << std::setw(5) << atoms[i].atom_index; // 7-11
        pdb_stream << " "; // 12
        if (atoms[i].atom.size() == 4) {
            pdb_stream << std::setw(4) << std::left << atoms[i].atom; // 13-16
        } else {
            pdb_stream << " ";
            pdb_stream << std::setw(3) << std::left << atoms[i].atom; // 13-16
        }
        pdb_stream << " "; // 17
        pdb_stream << std::setw(3) << std::right << atoms[i].residue; // 18-20
        pdb_stream << " "; // 21
        pdb_stream << atoms[i].chain; // 22
        pdb_stream << std::setw(4) << atoms[i].residue_index; // 23-26
        pdb_stream << "    "; // 27-30
        char buffer[16];
        fast_ftoa<1000, 3>(atoms[i].coordinate.x, buffer);
        pdb_stream << std::setw(8) << buffer; // 31-38
        fast_ftoa<1000, 3>(atoms[i].coordinate.y, buffer);
        pdb_stream << std::setw(8) << buffer; // 39-46
        fast_ftoa<1000, 3>(atoms[i].coordinate.z, buffer);
        pdb_stream << std::setw(8) << buffer; // 47-54
        pdb_stream << "  1.00"; // 55-60
        fast_ftoa<100, 2>(atoms[i].tempFactor, buffer);
        pdb_stream << std::setw(6) << buffer; // 61-66
        pdb_stream << "          "; // 67-76
        // First one character from atom
        pdb_stream << std::setw(2) << atoms[i].atom[0]; // 77-78
        pdb_stream << "  \n"; // 79-80
        if (i == (total-1)) {
            // TER
            // 1-6 Record name "TER   "
            // 7-11 Atom serial number.
            // 18-20 Residue name.
            // 22 Chain identifier.
            // 23-26 Residue sequence number.
            pdb_stream << "TER   " << std::setw(5) << atoms[i].atom_index + 1 << "      ";
            pdb_stream << std::setw(3) << std::right << atoms[i].residue;
            pdb_stream << " " << atoms[i].chain;
            pdb_stream << std::setw(4) << atoms[i].residue_index << std::endl;
        }
    }
}

int writeAtomCoordinatesToPDBFile(
    std::vector<AtomCoordinate>& atoms, std::string title, std::string pdb_path
) {
    std::ofstream pdb_file(pdb_path);
    if (!pdb_file) {
        return 1;
    }
    writeAtomCoordinatesToPDB(atoms, title, pdb_file);
    return 0;
}

std::vector< std::vector<AtomCoordinate> > splitAtomByResidue(
    const tcb::span<AtomCoordinate>& atomCoordinates
) {
    std::vector< std::vector<AtomCoordinate> > output;
    std::vector<AtomCoordinate> currentResidue;

    for (size_t i = 0; i < atomCoordinates.size(); i++) {
        if (i == 0) {
            currentResidue.push_back(atomCoordinates[i]);
        } else if (i != (atomCoordinates.size() - 1)) {
            if (atomCoordinates[i].residue_index == atomCoordinates[i-1].residue_index) {
                currentResidue.push_back(atomCoordinates[i]);
            } else {
                output.push_back(currentResidue);
                currentResidue.clear();
                currentResidue.push_back(atomCoordinates[i]);
            }
        } else {
            currentResidue.push_back(atomCoordinates[i]);
            output.push_back(currentResidue);
        }
    }

    return output;
}

std::vector<std::string> getResidueNameVector(
    const tcb::span<AtomCoordinate>& atomCoordinates
) {
    std::vector<std::string> output;
    // Unique residue names
    for (size_t i = 0; i < atomCoordinates.size(); i++) {
        if (i == 0) {
            output.push_back(atomCoordinates[i].residue);
        } else {
            if (atomCoordinates[i].residue_index != atomCoordinates[i-1].residue_index) {
                output.push_back(atomCoordinates[i].residue);
            }
        }
    }
    return output;
}

AtomCoordinate findFirstAtom(const std::vector<AtomCoordinate>& atoms, std::string atom_name) {
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.atom == atom_name) {
            return curr_atm;
        }
    }
    return AtomCoordinate();
}

void setAtomIndexSequentially(std::vector<AtomCoordinate>& atoms, int start) {
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms[i].atom_index = start + i;
    }
}

void removeAlternativePosition(std::vector<AtomCoordinate>& atoms) {
    // If there is an alternative position, remove it
    for (size_t i = 1; i < atoms.size(); i++) {
        if (atoms[i].atom == atoms[i-1].atom) {
            atoms.erase(atoms.begin() + i);
            i--;
        }
    }
}

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    std::vector<AtomCoordinate>& atoms, int residue_index
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index == residue_index) {
            output.emplace_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> getAtomsWithResidueIndiceRange(
    std::vector<AtomCoordinate>& atoms, int start, int end
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index >= start && curr_atm.residue_index < end) {
            output.emplace_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, int residue_index,
    std::vector<std::string> atomNames
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index == residue_index) {
            for (const std::string& atom_name : atomNames) {
                if (curr_atm.atom == atom_name) {
                    output.emplace_back(curr_atm);
                }
            }
        }
    }
    return output;
}

std::vector< std::vector<AtomCoordinate> > getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, std::vector<int> residue_index,
    std::vector<std::string> atomNames
) {
    std::vector<std::vector<AtomCoordinate>> output;
    for (int curr_index : residue_index) {
        output.emplace_back(getAtomsWithResidueIndex(atoms, curr_index, atomNames));
    }
    return output;
}

float RMSD(std::vector<AtomCoordinate>& atoms1, std::vector<AtomCoordinate>& atoms2) {
    // RMSD: Root Mean Square Deviation
    float sum = 0;
    // Sum of square of distance
    for (size_t i = 0; i < atoms1.size(); i++) {
        sum += pow(atoms1[i].coordinate.x - atoms2[i].coordinate.x, 2);
        sum += pow(atoms1[i].coordinate.y - atoms2[i].coordinate.y, 2);
        sum += pow(atoms1[i].coordinate.z - atoms2[i].coordinate.z, 2);
    }
    return sqrt(sum / atoms1.size());
}

std::vector<AtomCoordinate> _subsetAtomVectorWithIndices(
    std::vector<AtomCoordinate>& atoms,
    std::pair<size_t, size_t>& indices
) {
    std::vector<AtomCoordinate> output;
    for (size_t i = indices.first; i < indices.second; i++) {
        output.push_back(atoms[i]);
    }
    return output;
}

void _splitAtomVectorWithIndices(
    std::vector<AtomCoordinate>& atoms,
    std::vector< std::pair<size_t, size_t> >& indices,
    std::vector< std::vector<AtomCoordinate> >& output
) {
    if (indices.size() == 0) {
        output.push_back(atoms);
    } else {
        for (size_t i = 0; i < indices.size(); i++) {
            output.push_back(_subsetAtomVectorWithIndices(atoms, indices[i]));
        }
    }
}

/**
 * @brief Identify discontinuous regions in atom coordinate vector and return
 *        vector of indices of the start and end of each region.
 *        start: inclusive, end: exclusive [start, end)
 * @param atoms
 * @param mode
 * @return std::vector<std::pair<size_t, size_t>>
 */
std::vector< std::pair<size_t, size_t> > identifyChains(const std::vector<AtomCoordinate>& atoms) {
    std::vector< std::pair<size_t, size_t> > output;
    size_t start = 0;
    // Split by chain
    for (size_t i = 1; i < atoms.size(); i++) {
        if (atoms[i].chain != atoms[i - 1].chain) {
            // Ensure that the new fragment starts with "N"
            if (atoms[i].atom == "N") {
                output.emplace_back(start, i);
                start = i;
            } else {
                // Find the first "N" atom
                for (size_t j = i; j < atoms.size(); j++) {
                    if (atoms[j].atom == "N") {
                        // Ignore fragment between i and j
                        output.emplace_back(start, i);
                        start = j;
                        break;
                    }
                }
                // Set i to j
                i = start;
            }
        }
    }
    // Add the last fragment
    output.emplace_back(start, atoms.size());
    return output;

}

/**
 * @brief Identify discontinuous residue indices in atom coordinate vector and return
 *        coordinates that have the same chain
 * @param atoms
 * @return std::vector< std::pair<size_t, size_t> >
 */
std::vector<std::pair<size_t, size_t>> identifyDiscontinousResInd(
    const std::vector<AtomCoordinate>& atoms,
    size_t chain_start,
    size_t chain_end
) {
    std::vector<std::pair<size_t, size_t>> output;
    // Extract N atoms only within chain
    std::vector<std::pair<size_t, int>> N_indices;
    for (size_t i = chain_start; i < chain_end; i++) {
        if (atoms[i].atom == "N") {
            N_indices.emplace_back(i, atoms[i].residue_index);
        }
    }
    // Identify discontinuous regions
    size_t start = N_indices[0].first;
    for (size_t i = 1; i < N_indices.size(); i++) {
        if (N_indices[i].second - N_indices[i - 1].second > 1) {
            output.emplace_back(start, N_indices[i].first);
            start = N_indices[i].first;
        }
    }
    // Add the last fragment
    output.emplace_back(start, chain_end);
    return output;
}
