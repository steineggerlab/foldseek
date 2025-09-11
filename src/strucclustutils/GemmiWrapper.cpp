//
// Created by Martin Steinegger on 6/7/21.
//
#include "GemmiWrapper.h"
#include "mmread.hpp"
#ifdef HAVE_ZLIB
#include "gz.hpp"
#endif
#include "input.hpp"
#include "foldcomp.h"
#include "cif.hpp"

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>

#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>

size_t decompressZstdBuffer(const void* src, size_t srcSize, char*& dst) {
    if (!src || srcSize == 0) {
        return 0;
    }

    const size_t frameSize = ZSTD_getFrameContentSize(src, srcSize);
    if (frameSize != ZSTD_CONTENTSIZE_ERROR && frameSize != ZSTD_CONTENTSIZE_UNKNOWN) {
        dst = new char[frameSize];
        const size_t ret = ZSTD_decompress(dst, frameSize, src, srcSize);
        if (ZSTD_isError(ret)) {
            delete[] dst;
            dst = nullptr;
            return 0;
        }
        return ret;
    }

    ZSTD_DStream* dctx = ZSTD_createDStream();
    if (!dctx) {
        return 0;
    }

    const size_t inChunk = ZSTD_DStreamInSize();
    const size_t outChunk = ZSTD_DStreamOutSize();

    std::vector<char> out;
    out.reserve(outChunk * 4);

    const char* ip = static_cast<const char*>(src);
    size_t remaining = srcSize;
    while (remaining) {
        const size_t toRead = remaining < inChunk ? remaining : inChunk;
        ZSTD_inBuffer inBuf { const_cast<char*>(ip), toRead, 0 };
        remaining -= toRead;
        ip += toRead;

        while (inBuf.pos < inBuf.size)  {
            const size_t oldSize = out.size();
            out.resize(oldSize + outChunk);

            ZSTD_outBuffer outBuf { out.data() + oldSize, outChunk, 0 };
            const size_t ret = ZSTD_decompressStream(dctx, &outBuf, &inBuf);
            if (ZSTD_isError(ret)) {
                ZSTD_freeDStream(dctx);
                return 0;
            }
            out.resize(oldSize + outBuf.pos);
        }
    }
    ZSTD_freeDStream(dctx);

    dst = new char[out.size()];
    memcpy(dst, out.data(), out.size());
    return out.size();
}

size_t decompressZstdFile(const std::string& path, char*& dst) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0 || st.st_size <= 0) {
        return 0;
    }

    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) {
        return 0;
    }

    const size_t fsize = static_cast<size_t>(st.st_size);
    std::vector<char> compressed(fsize);
    if (fread(compressed.data(), 1, fsize, fp) != fsize) {
        fclose(fp);
        return 0;
    }
    fclose(fp);

    return decompressZstdBuffer(compressed.data(), fsize, dst);
}

GemmiWrapper::GemmiWrapper(){
    fixupBuffer = NULL;
}

// adapted from Gemmi find_tabulated_residue
char threeToOneAA(const std::string& three) {
    if (three.size() != 3) {
        return 'X';
    }
#define ID(s) (s[0] << 16 | s[1] << 8 | s[2])
    switch (ID(three.c_str())) {
        case ID("ALA"): return 'A';
        case ID("ARG"): return 'R';
        case ID("ASN"): return 'N';
        case ID("ABA"): return 'A';
        case ID("ASP"): return 'D';
        case ID("ASX"): return 'B';
        case ID("CYS"): return 'C';  // also BUF
        case ID("CSH"): return 'S';
        case ID("GLN"): return 'Q';
        case ID("GLU"): return 'E';
        case ID("GLX"): return 'Z';
        case ID("GLY"): return 'G';  // also BUF
        case ID("HIS"): return 'H';
        case ID("ILE"): return 'I';
        case ID("LEU"): return 'L';
        case ID("LYS"): return 'K';
        case ID("MET"): return 'M';
        case ID("MSE"): return 'M';
        case ID("ORN"): return 'A';
        case ID("PHE"): return 'F';
        case ID("PRO"): return 'P';
        case ID("SER"): return 'S';
        case ID("THR"): return 'T';
        case ID("TRY"): // fall-through - synonym for TRP
        case ID("TRP"): return 'W';
        case ID("TYR"): return 'Y';
        case ID("UNK"): return 'X';
        case ID("VAL"): return 'V';
        case ID("SEC"): return 'C'; // mapped to U in gemmi, C in previous FS
        case ID("PYL"): return 'O';
        case ID("SEP"): return 'S';
        case ID("TPO"): return 'T';
        case ID("PCA"): return 'E';
        case ID("CSO"): return 'C';
        case ID("PTR"): return 'Y';
        case ID("KCX"): return 'K';
        case ID("CSD"): return 'C';
        case ID("LLP"): return 'K';
        case ID("CME"): return 'C';
        case ID("MLY"): return 'K';
        case ID("DAL"): return 'A';
        case ID("TYS"): return 'Y';
        case ID("OCS"): return 'C';
        case ID("M3L"): return 'K';
        case ID("FME"): return 'M';
        case ID("ALY"): return 'K';
        case ID("HYP"): return 'P';
        case ID("CAS"): return 'C';
        case ID("CRO"): return 'T';
        case ID("CSX"): return 'C';
        case ID("DPR"): return 'P';  // also BUF
        case ID("DGL"): return 'E';
        case ID("DVA"): return 'V';
        case ID("CSS"): return 'C';
        case ID("DPN"): return 'F';
        case ID("DSN"): return 'S';
        case ID("DLE"): return 'L';
        case ID("HIC"): return 'H';
        case ID("NLE"): return 'L';
        case ID("MVA"): return 'V';
        case ID("MLZ"): return 'K';
        case ID("CR2"): return 'G';
        case ID("SAR"): return 'G';
        case ID("DAR"): return 'R';
        case ID("DLY"): return 'K';
        case ID("YCM"): return 'C';
        case ID("NRQ"): return 'M';
        case ID("CGU"): return 'E';
        case ID("0TD"): return 'D';
        case ID("MLE"): return 'L';
        case ID("DAS"): return 'D';
        case ID("DTR"): return 'W';
        case ID("CXM"): return 'M';
        case ID("TPQ"): return 'Y';
        case ID("DCY"): return 'C';
        case ID("DSG"): return 'N';
        case ID("DTY"): return 'Y';
        case ID("DHI"): return 'H';
        case ID("MEN"): return 'N';
        case ID("DTH"): return 'T';
        case ID("SAC"): return 'S';
        case ID("DGN"): return 'Q';
        case ID("AIB"): return 'A';
        case ID("SMC"): return 'C';
        case ID("IAS"): return 'D';
        case ID("CIR"): return 'R';
        case ID("BMT"): return 'T';
        case ID("DIL"): return 'I';
        case ID("FGA"): return 'E';
        case ID("PHI"): return 'F';
        case ID("CRQ"): return 'Q';
        case ID("SME"): return 'M';
        case ID("GHP"): return 'G';
        case ID("MHO"): return 'M';
        case ID("NEP"): return 'H';
        case ID("TRQ"): return 'W';
        case ID("TOX"): return 'W';
        case ID("ALC"): return 'A';
        //   case ID("3FG"): return ' ';
        case ID("SCH"): return 'C';
        case ID("MDO"): return 'A';
        case ID("MAA"): return 'A';
        case ID("GYS"): return 'S';
        case ID("MK8"): return 'L';
        case ID("CR8"): return 'H';
        case ID("KPI"): return 'K';
        case ID("SCY"): return 'C';
        case ID("DHA"): return 'S';
        case ID("OMY"): return 'Y';
        case ID("CAF"): return 'C';
        case ID("0AF"): return 'W';
        case ID("SNN"): return 'N';
        case ID("MHS"): return 'H';
        //   case ID("MLU"): return ' ';
        case ID("SNC"): return 'C';
        case ID("PHD"): return 'D';
        case ID("B3E"): return 'E';
        case ID("MEA"): return 'F';
        case ID("MED"): return 'M';
        case ID("OAS"): return 'S';
        case ID("GL3"): return 'G';
        case ID("FVA"): return 'V';
        case ID("PHL"): return 'F';
        case ID("CRF"): return 'T';
        //   case ID("OMZ"): return ' ';
        case ID("BFD"): return 'D';
        case ID("MEQ"): return 'Q';
        case ID("DAB"): return 'A';
        case ID("AGM"): return 'R';
        // added from previous FS code
        case ID("4BF"): return 'Y';
        case ID("B3A"): return 'A';
        case ID("B3D"): return 'D';
        case ID("B3K"): return 'K';
        case ID("B3Y"): return 'Y';
        case ID("BAL"): return 'A';
        case ID("DBZ"): return 'A';
        case ID("GPL"): return 'K';
        case ID("HSK"): return 'H';
        case ID("HY3"): return 'P';
        case ID("HZP"): return 'P';
        case ID("KYN"): return 'W';
        case ID("MGN"): return 'Q';
    }
    return 'X';
#undef ID
}

std::unordered_map<std::string, int> getEntityTaxIDMapping(gemmi::cif::Document& doc) {
    std::unordered_map<std::string, int> entity_to_taxid;
    static const std::vector<std::pair<std::string, std::string>> loops_with_taxids = {
        { "_entity_src_nat.", "?pdbx_ncbi_taxonomy_id"},
        { "_entity_src_gen.", "?pdbx_gene_src_ncbi_taxonomy_id"},
        { "_pdbx_entity_src_syn.", "?ncbi_taxonomy_id"}
    };
    for (gemmi::cif::Block& block : doc.blocks) {
        for (auto&& [loop, taxid] : loops_with_taxids) {
            for (auto row : block.find(loop, {"entity_id", taxid})) {
                if (row.has2(1) == false) {
                    continue;
                }
                std::string entity_id = gemmi::cif::as_string(row[0]);
                if (entity_to_taxid.find(entity_id) != entity_to_taxid.end()) {
                    continue;
                }
                const char* endptr = NULL;
                int taxId = gemmi::no_sign_atoi(row[1].c_str(), &endptr);
                if (endptr != NULL && *endptr == '\0') {
                    entity_to_taxid.emplace(entity_id, taxId);
                }
            }
        }
    }
    return entity_to_taxid;
}

/*std::unordered_map<std::string, unsigned int> getResidueIndices(gemmi::cif::Document & doc) {
    std::unordered_map<std::string, unsigned int> residueIndices;
    for (gemmi::cif::Block & block : doc.blocks) {
        for (auto row : block.find("_atom_site.", {"label_asym_id",
                                                   "?label_seq_id", 
                                                   "?auth_seq_id",
                                                   "?auth_asym_id"})) {
 
            std::string resLabel;
 

            if (row.has2(4) == true) {
 

                //If auth chain label exists, start map key with it
 

                resLabel = gemmi::cif::as_string(row[4]);
 

            }
 

            else {
 

                //If auth chain label does not exist, start map key with asym label
 

                resLabel = gemmi::cif:as_string(row[0]);
 

            }
 

            if (row.has2(3) == true) {
 

                //If auth seq id exists, append to chain id as the map key
 

                resLabel.append(gemmi::cif::as_string(row[3]));
 

            }
 

            else {
 

                //If auth seq id does not exist, append label seq id
 

                resLabel.append(gemmi::cif::as_string(row[2]));
 

            }
 

            if (residueIndices.find(resLabel) != residueIndices.end()) {
 

                //skip residue label if already mapped
 

                continue;
 

            }
 

            const char* endptr = NULL;
 

            unsigned int resId = gemmi::no_sign_atoi(row[1].c_str(), &endptr);
 

            if (endptr != NULL && *endptr == '\0') {
 

                residueIndices.emplace(resLabel, resId);
 

            }
 

        }
 

    }
 

    return residueIndices;
 

}*/
 

GemmiWrapper::Format mapFormat(gemmi::CoorFormat format) {
    switch (format) {
        case gemmi::CoorFormat::Pdb:
            return GemmiWrapper::Format::Pdb;
        case gemmi::CoorFormat::Mmcif:
            return GemmiWrapper::Format::Mmcif;
        case gemmi::CoorFormat::Mmjson:
            return GemmiWrapper::Format::Mmjson;
        case gemmi::CoorFormat::ChemComp:
            return GemmiWrapper::Format::ChemComp;
        default:
            return GemmiWrapper::Format::Unknown;
    }
}

bool GemmiWrapper::load(const std::string& filename, bool saveResIndex, Format format) {
    if ((format == Format::Foldcomp) || (format == Format::Detect && gemmi::iends_with(filename, ".fcz"))) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            return false;
        }
        return loadFoldcompStructure(in, filename);
    }
    try {
        if (gemmi::iends_with(filename, ".zstd") || gemmi::iends_with(filename, ".zst")) {
            std::string name = gemmi::path_basename(filename, { ".zstd", ".zst" });
            char* out;
            size_t len = decompressZstdFile(filename, out);
            if (format == Format::Detect) {
                format = mapFormat(gemmi::coor_format_from_ext(name));
            }
            if (format == Format::Unknown) {
                format = Format::Pdb;
            }
            return loadFromBuffer(out, len, name, saveResIndex, format);
        }

#ifdef HAVE_ZLIB
        gemmi::MaybeGzipped infile(filename);
#else
        gemmi::BasicInput infile(filename);
#endif
        if (format == Format::Detect) {
            format = mapFormat(gemmi::coor_format_from_ext(infile.basepath()));
        }
        gemmi::Structure st;
        std::unordered_map<std::string, int> entity_to_tax_id;
        //std::unordered_map<std::string, unsigned int> residueIndices;
        switch (format) {
            case Format::Mmcif: {
                gemmi::CharArray mem = gemmi::read_into_buffer(infile);
                char* data = mem.data();
                size_t dataSize = mem.size();

                // hack to fix broken _citation.title in AF3
                const char target0[] = "Accurate structure prediction of biomolecular interactions with AlphaFold 3\n";
                size_t target0Len = sizeof(target0) - 1;
                const char target1[] = "_citation.title";
                size_t target1Len = sizeof(target1) - 1;
                char* it = std::search(data, data + dataSize, target0, target0 + target0Len);
                if (it != data + dataSize) {
                    while (it > data && *(it - 1) != '\n') {
                        it--;
                    }
                    if (strncmp(it, target1, target1Len) == 0) {
                        it[0] = '#';
                        it[1] = ' ';
                    }
                }

                gemmi::cif::Document doc = gemmi::cif::read_memory(mem.data(), mem.size(), infile.path().c_str());
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                //residueIndices = getResidueIndices(doc);
                st = gemmi::make_structure(doc);
                break;
            }
            case Format::Mmjson: {
                gemmi::cif::Document doc = gemmi::cif::read_mmjson(infile);
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure(doc);
                //residueIndices = getResidueIndices(doc);
                break;
            }
            case Format::ChemComp: {
                gemmi::cif::Document doc = gemmi::cif::read(infile);
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure_from_chemcomp_doc(doc);
                //residueIndices = getResidueIndices(doc);
                break;
            }
            default:
                st = gemmi::read_pdb(infile);
        }
        updateStructure((void*) &st, filename, entity_to_tax_id, saveResIndex);
    } catch (...) {
        return false;
    }
    return true;
}

// https://stackoverflow.com/questions/1448467/initializing-a-c-stdistringstream-from-an-in-memory-buffer/1449527
struct OneShotReadBuf : public std::streambuf
{
    OneShotReadBuf(char* s, std::size_t n)
    {
        setg(s, s, s + n);
    }
};
bool GemmiWrapper::loadFromBuffer(const char * buffer, size_t bufferSize, const std::string& name, bool saveResIndex, GemmiWrapper::Format format) {
    if ((format == Format::Foldcomp) || (format == Format::Detect && (bufferSize > MAGICNUMBER_LENGTH && strncmp(buffer, MAGICNUMBER, MAGICNUMBER_LENGTH) == 0))) {
        OneShotReadBuf buf((char *) buffer, bufferSize);
        std::istream istr(&buf);
        if (!istr) {
            return false;
        }
        return loadFoldcompStructure(istr, name);
    }
    try {
        std::string newName = name;
        const char* newBuffer = buffer;
        size_t newBufferSize = bufferSize;
        bool zstd = false;
        if (gemmi::iends_with(name, ".zstd") || gemmi::iends_with(name, ".zst")) {
            newName = gemmi::path_basename(name, { ".zstd", ".zst" });
            char* tmpBuffer;
            newBufferSize = decompressZstdBuffer(buffer, bufferSize, tmpBuffer);
            newBuffer = tmpBuffer;
        }

#ifdef HAVE_ZLIB
        gemmi::MaybeGzipped infile(newName);
#else
        gemmi::BasicInput infile(newName);
#endif
        if (format == Format::Detect) {
            format = mapFormat(gemmi::coor_format_from_ext(infile.basepath()));
        }

        gemmi::Structure st;
        std::unordered_map<std::string, int> entity_to_tax_id;
        //std::unordered_map<std::string, unsigned int> residueIndices;
        switch (format) {
            case Format::Pdb:
                st = gemmi::pdb_impl::read_pdb_from_stream(gemmi::MemoryStream(newBuffer, newBufferSize), newName, gemmi::PdbReadOptions());
                break;
            case Format::Mmcif: {
                const char* targetBuffer = newBuffer;
                // hack to fix broken _citation.title in AF3
                const char target0[] = "Accurate structure prediction of biomolecular interactions with AlphaFold 3\n";
                size_t target0Len = sizeof(target0) - 1;
                const char target1[] = "_citation.title";
                size_t target1Len = sizeof(target1) - 1;
                const char* it = std::search(targetBuffer, targetBuffer + newBufferSize, target0, target0 + target1Len);
                if (it != targetBuffer + newBufferSize) {
                    if (fixupBuffer == NULL) {
                        fixupBufferSize = newBufferSize;
                        fixupBuffer = (char*)malloc(fixupBufferSize);
                    } else if (newBufferSize > fixupBufferSize) {
                        fixupBufferSize = bufferSize * 1.5;
                        fixupBuffer = (char*)realloc(fixupBuffer, fixupBufferSize);
                    }
                    memcpy(fixupBuffer, targetBuffer, bufferSize);
                    while (it > targetBuffer && *(it - 1) != '\n') {
                        it--;
                    }
                    if (strncmp(it, target1, target1Len) == 0) {
                        *(fixupBuffer + (it - targetBuffer)) = '#';
                        *(fixupBuffer + (it - targetBuffer) + 1) = ' ';
                    }
                    targetBuffer = fixupBuffer;
                }
                gemmi::cif::Document doc = gemmi::cif::read_memory(targetBuffer, newBufferSize, newName.c_str());
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                //residueIndices = getResidueIndices(doc);
                st = gemmi::make_structure(doc);
                break;
            }
            case Format::Mmjson: {
                char* bufferCopy = (char*)malloc(newBufferSize + 1 * sizeof(char));
                if (bufferCopy == NULL) {
                    return false;
                }
                if (memcpy(bufferCopy, newBuffer, newBufferSize) == NULL) {
                    free(bufferCopy);
                    return false;
                }
                bufferCopy[newBufferSize] = '\0';
                gemmi::cif::Document doc = gemmi::cif::read_mmjson_insitu(bufferCopy, newBufferSize, newName);
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                //residueIndices = getResidueIndices(doc);
                st = gemmi::make_structure(doc);
                free(bufferCopy);
                break;
            }
            case Format::ChemComp: {
                gemmi::cif::Document doc = gemmi::cif::read_memory(newBuffer, newBufferSize, newName.c_str());
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                //residueIndices = getResidueIndices(doc);
                st = gemmi::make_structure_from_chemcomp_doc(doc);
                break;
            }
            default:
                return false;
        }
        updateStructure((void*) &st, newName, entity_to_tax_id, saveResIndex);
    } catch (...) {
        return false;
    }
    return true;
}

bool GemmiWrapper::loadFoldcompStructure(std::istream& stream, const std::string& filename) {
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
    seq3di.clear();
    taxIds.clear();
    title.clear();
    chain.clear();
    names.clear();
    chainNames.clear();
    modelIndices.clear();
    ca.clear();
    ca_bfactor.clear();
    c.clear();
    cb.clear();
    n.clear();
    ami.clear();
    title.append(fc.strTitle);
    names.push_back(filename);
    const AtomCoordinate& first = coordinates[0];
    chainNames.push_back(first.chain);
    modelCount = 1;
    modelIndices.push_back(modelCount);
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
            ami.push_back(threeToOneAA(atom.residue));
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

void GemmiWrapper::updateStructure(void * void_st, const std::string& filename, std::unordered_map<std::string, int>& entity_to_tax_id, bool saveResIndex) {
    gemmi::Structure * st = (gemmi::Structure *) void_st;

    title.clear();
    chain.clear();
    names.clear();
    chainNames.clear();
    modelIndices.clear();
    modelCount = 0;
    ca.clear();
    ca_bfactor.clear();
    c.clear();
    cb.clear();
    n.clear();
    ami.clear();
    taxIds.clear();
    resIds.clear();
    title.append(st->get_info("_struct.title"));
    size_t currPos = 0;
    for (gemmi::Model& model : st->models){
        modelCount++;
        for (gemmi::Chain& ch : model.chains) {
            size_t chainStartPos = currPos;
            size_t pos = filename.find_last_of("\\/");
            std::string name = (std::string::npos == pos)
                               ? filename
                               : filename.substr(pos+1, filename.length());
            //name.push_back('_');
            chainNames.push_back(ch.name);
            char* rest;
            errno = 0;
            unsigned int modelNumber = strtoul(model.name.c_str(), &rest, 10);
            if ((rest != model.name.c_str() && *rest != '\0') || errno == ERANGE) {
                modelIndices.push_back(modelCount);
            }else{
                modelIndices.push_back(modelNumber);
            }

            names.push_back(name);
            int taxId = -1;
            for (const gemmi::Residue &res : ch.first_conformer()) {
                // only consider polymers and unknowns
                if (res.entity_type != gemmi::EntityType::Polymer && res.entity_type != gemmi::EntityType::Unknown) {
                    continue;
                }
                // only consider ATOM and HETATM
                if (res.het_flag == '\0') {
                    continue;
                }

                Vec3 ca_atom = {NAN, NAN, NAN};
                Vec3 cb_atom = {NAN, NAN, NAN};
                Vec3 n_atom  = {NAN, NAN, NAN};
                Vec3 c_atom  = {NAN, NAN, NAN};
                float ca_atom_bfactor = 0.0f;
                bool hasCA = false;
                for (const gemmi::Atom &atom : res.atoms) {
                    if (atom.name == "CA") {
                        ca_atom.x = atom.pos.x;
                        ca_atom.y = atom.pos.y;
                        ca_atom.z = atom.pos.z;
                        ca_atom_bfactor = atom.b_iso;
                        hasCA = true;
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
                if (hasCA == false) {
                    continue;
                }
                ca_bfactor.push_back(ca_atom_bfactor);
                ca.push_back(ca_atom);
                cb.push_back(cb_atom);
                n.push_back(n_atom);
                c.push_back(c_atom);

                ami.push_back(threeToOneAA(res.name));

                if (saveResIndex) {
                    if (!res.label_seq){
                        int tmpResId = static_cast<int>(res.seqid.num);
                        unsigned int resId = unsigned(tmpResId);
                        resIds.push_back(resId);
                    }
                    else {
                        int tmpResId = static_cast<int>(res.label_seq);
                        unsigned int resId = unsigned(tmpResId);
                        resIds.push_back(resId);
                    }
                    //std::cout << ch.name << " " << resId << std::endl;
                }

                if (taxId == -1) {
                    auto it = entity_to_tax_id.find(res.entity_id);
                    if (it != entity_to_tax_id.end()) {
                        taxId = it->second;
                    }
                }
                currPos++;
            }
            taxIds.push_back(taxId == -1 ? 0 : taxId);
            chain.push_back(std::make_pair(chainStartPos, currPos));
        }
    }
}
