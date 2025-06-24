/**
 * File: structure_reader.cpp
 * Project: foldcomp
 * Created: 2022-08-04 15:02:09
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "src".
 * ---
 * Last Modified: 2022-10-18 17:50:57
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2022 Hyunbin Kim, All rights reserved
 */
#include "structure_reader.h"

#include <stdexcept>
#include <zlib.h>
#include <iostream>
// Gemmi
#include "gemmi/gz.hpp"
#include "gemmi/input.hpp"
#include "gemmi/mmread.hpp"

/**
 * @brief StructureReader::updateStructure
 * Reads a structure and updates the AtomCoordinate vector.
 * Inner function of load() and loadFromBuffer().
 * @param void_st
 * @param filename
 */
void StructureReader::updateStructure(void* void_st, const std::string& filename) {
    gemmi::Structure* st = (gemmi::Structure* ) void_st;

    this->title.clear();
    this->atoms.clear();
    // cif
    this->title = st->get_info("_entry.id");
    // pdb
    if (this->title.empty()) {
        this->title = st->get_info("_struct.title");
    }
    // else
    if (this->title.empty()) {
        this->title = filename;
    }

    for (gemmi::Model& model : st->models) {
        for (gemmi::Chain& ch : model.chains) {
            for (gemmi::Residue& res : ch.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    AtomCoordinate ac = AtomCoordinate(
                        atom.name, res.name, ch.name, atom.serial, (int)res.seqid.num,
                        (float)atom.pos.x, (float)atom.pos.y, (float)atom.pos.z
                    );
                    ac.tempFactor = atom.b_iso;
                    this->atoms.push_back(ac);
                }
            }
        }
    }
}

/**
 * @brief StructureReader::loadFromBuffer
 * Reads a structure from a buffer and updates the AtomCoordinate vector.
 * Originally from GemmiWrapper.cpp in github.com/steineggerlab/foldseek
 * Original author: Martin Steinegger
 * @param buffer
 * @param bufferSize
 * @param name
 * @return true
 * @return false
 */
bool StructureReader::loadFromBuffer(const char* buffer, size_t bufferSize, const std::string& name) {
    try {
        gemmi::MaybeGzipped infile(name);
        gemmi::Structure st;
        // If infile is a compressed file, we need to uncompress buffer using zlib
        if (infile.is_compressed()) {
            char* uncompBuffer;
            size_t uncompBufferSize;
            int ret = uncompressBuffer(const_cast<const char**>(&uncompBuffer), &uncompBufferSize, buffer, bufferSize);
            if (ret != 0) {
                return false;
            }
            st = gemmi::read_structure_from_char_array(uncompBuffer, uncompBufferSize, name);
            free(uncompBuffer);
        } else {
            st = gemmi::read_structure_from_char_array(const_cast<char*>(buffer), bufferSize, name);
        }

        updateStructure((void*)&st, name);
    }
    catch (std::runtime_error& e) {
        return false;
    }
    return true;
}

/**
 * @brief StructureReader::openStructure
 * Originally from GemmiWrapper.cpp in github.com/steineggerlab/foldseek
 * Original author: Martin Steinegger
 * @param filename
 * @return gemmi::Structure
 */
gemmi::Structure openStructure(const std::string& filename) {
    gemmi::MaybeGzipped infile(filename);
    gemmi::CoorFormat format = gemmi::coor_format_from_ext(infile.basepath());
    if (format != gemmi::CoorFormat::Unknown && format != gemmi::CoorFormat::Unknown) {
        return gemmi::read_structure(infile, format);
    }
    else {
        return gemmi::read_structure(infile, gemmi::CoorFormat::Pdb);
    }
}

/**
 * @brief StructureReader::load
 * Originally from GemmiWrapper.cpp in github.com/steineggerlab/foldseek
 * Original author: Martin Steinegger
 * @param filename
 * @return true
 * @return false
 */
bool StructureReader::load(const std::string& filename){
    try {
        gemmi::Structure st = openStructure(filename);
        updateStructure((void*)&st, filename);
    }
    catch (std::runtime_error& e) {
        return false;
    }
    return true;
}

bool StructureReader::readBackboneAtoms(std::vector<AtomCoordinate>& backboneAtoms) {
    backboneAtoms.clear();
    for (AtomCoordinate& ac : this->atoms) {
        if (ac.atom == "CA" || ac.atom == "C" || ac.atom == "N") {
            backboneAtoms.push_back(ac);
        }
    }
    return true;
}

bool StructureReader::readAllAtoms(std::vector<AtomCoordinate>& allAtoms){
    allAtoms.clear();
    // copy all atoms
    for (AtomCoordinate& ac : this->atoms) {
        allAtoms.push_back(ac);
    }
    return true;
}

int uncompressBuffer(
    const char** uncompBuffer, size_t* uncompBufferSize,
    const char* origBuffer, size_t origBufferSize
) {
    const unsigned int CHUNK = 128 * 1024;
    unsigned char out[CHUNK];
    z_stream strm;
    memset(&strm, 0, sizeof(z_stream));
    int status = inflateInit2(&strm, 15 | 32);
    if (status < 0) {
        return -1;
    }
    strm.avail_in = origBufferSize;
    strm.next_in = (unsigned char*)origBuffer;
    size_t bufferSize = 1024 * 1024;
    char* buffer = (char*)malloc(bufferSize);
    char* current = buffer;
    do {
        unsigned have;
        strm.avail_out = CHUNK;
        strm.next_out = out;
        int err = inflate(&strm, Z_NO_FLUSH);
        switch (err) {
        case Z_OK:
        case Z_STREAM_END:
        case Z_BUF_ERROR:
            break;
        default:
            inflateEnd(&strm);
            free(buffer);
            return -1;
        }
        have = CHUNK - strm.avail_out;
        if (current + have > buffer + bufferSize) {
            size_t offset = current - buffer;
            bufferSize *= 2;
            buffer = (char*)realloc(buffer, bufferSize);
            current = buffer + offset;
        }
        memcpy(current, out, have);
        current += have;

    } while (strm.avail_out == 0);
    inflateEnd(&strm);
    *uncompBuffer = buffer;
    *uncompBufferSize = current - buffer;
    return 0;
}