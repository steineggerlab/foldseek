#ifndef PATHHASHER_H
#define PATHHASHER_H

#include <string>
#include <cstdint>

class PathHasher {
public:
    // Encodes a 64-bit hash to base62 string (alphanumeric only)
    static std::string hashToBase62(uint64_t hash);

    // Hashes a full path string and returns base62 encoded hash
    static std::string hashPath(const std::string& path);

private:
    static const char BASE62_ALPHABET[];
    static const int BASE62_BASE = 62;
};

#endif // PATHHASHER_H
