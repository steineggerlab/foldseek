#ifndef STRUCTURE_UTIL_H
#define STRUCTURE_UTIL_H

#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "LocalParameters.h"
#include "BaseMatrix.h"
#include "structureto12st.h"

#include <vector>

class StructureUtil {
public:
    static bool is3Di12StDb(int dbtype) {
        return (DBReader<unsigned int>::getExtendedDbtype(dbtype)
                & LocalParameters::DBTYPE_EXTENDED_3DI_12ST) != 0;
    }

    static inline void split3Di12St(const char *src, size_t len,
                                     std::vector<char> &seq3di,
                                     std::vector<char> &seq12st,
                                     const BaseMatrix &subMat3Di,
                                     const BaseMatrix &subMat12St) {
        if (seq3di.size() < len) {
            seq3di.resize(len);
        }
        if (seq12st.size() < len) {
            seq12st.resize(len);
        }
        for (size_t i = 0; i < len; ++i) {
            unsigned char val = static_cast<unsigned char>(src[i]);
            unsigned char state3di = static_cast<unsigned char>(val / Alphabet12St::STATE_CNT);
            unsigned char state12st = static_cast<unsigned char>(val % Alphabet12St::STATE_CNT);
            seq3di[i] = subMat3Di.num2aa[state3di];
            seq12st[i] = subMat12St.num2aa[state12st];
        }
    }

    static std::string getIndexWithSuffix(std::string db, const std::string &suffix) {
        if (Util::endsWith(".idx", db)) {
            db = db.substr(0, db.length() - 4);
        } else {
            return db + suffix;
        }
        db.append(suffix);
        std::string index = PrefilteringIndexReader::searchForIndex(db);
        if (index != "") {
            return index;
        }
        return db;
    }
};

#endif
