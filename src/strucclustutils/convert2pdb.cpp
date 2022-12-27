#include <cstring>
#include <cstdio>

#include "LocalParameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Coordinate16.h"

const char* threeLetterLookup[26] = { "ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "XLE", "LYS", "LEU", "MET", "ASN", "PYL", "PRO", "GLN", "ARG", "SER", "THR", "SEC", "VAL", "TRP", "XAA", "TYR", "GLX" };
int convert2pdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> db(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> db_header(par.hdr1.c_str(), par.hdr1Index.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db_header.open(DBReader<unsigned int>::NOSORT);

    std::string dbCa = par.db1 + "_ca";
    std::string dbCaIndex = par.db1 + "_ca.index";
    DBReader<unsigned int> db_ca(dbCa.c_str(), dbCaIndex.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db_ca.open(DBReader<unsigned int>::NOSORT);

    FILE* handle = fopen(par.db2.c_str(), "w");
    if(handle == NULL) {
        perror(par.db2.c_str());
        EXIT(EXIT_FAILURE);
    }

    Coordinate16 coords;
    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";
    for(size_t i = 0; i < db.getSize(); i++){
        unsigned int key = db.getDbKey(i);
        unsigned int headerId = db_header.getId(key);
        const char* headerData = db_header.getData(headerId, 0);
        const size_t headerLen = db_header.getEntryLen(headerId) - 2;

        unsigned int seqId = db.getId(key);
        const char* seqData = db.getData(seqId, 0);
        const size_t seqLen = std::max(db.getEntryLen(seqId), (size_t)2) - 2;

        unsigned int caId = db_ca.getId(key);
        const char* caData = db_ca.getData(caId, 0);
        const size_t caLen = db_ca.getEntryLen(caId);

        float* ca = coords.read(caData, seqLen, caLen);

        fprintf(handle, "MODEL % 8d\n", key);
        int remainingHeader = headerLen;
        fprintf(handle, "TITLE     %.*s\n",  std::min(70, (int)remainingHeader), headerData);
        remainingHeader -= 70;
        int continuation = 2;
        while (remainingHeader > 0) {
            fprintf(handle, "TITLE  % 3d%.*s\n", continuation, std::min(70, (int)remainingHeader), headerData + (headerLen - remainingHeader));
            remainingHeader -= 70;
            continuation++;
        }
        for (size_t j = 0; j < seqLen; ++j) {
            // make AA upper case
            char aa = seqData[j] & ~0x20;
            if (aa == '\0') {
                aa = 'X';
            }
            const char* aa3 = threeLetterLookup[(int)(aa - 'A')];
            fprintf(handle, "ATOM  % 5d  CA  %s A% 4d% 12.3f% 8.3f% 8.3f\n", (int)(j + 1), aa3, int(j + 1), ca[j], ca[j + (1 * seqLen)], ca[j + (2 * seqLen)]);
        }
        fprintf(handle, "ENDMDL\n");
    }
    if (fclose(handle) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << par.db2 << "\n";
        EXIT(EXIT_FAILURE);
    }
    db_ca.close();
    db_header.close();
    db.close();

    return EXIT_SUCCESS;
}
