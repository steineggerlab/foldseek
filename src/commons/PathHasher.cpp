#include "PathHasher.h"
#include <functional>

const char PathHasher::BASE62_ALPHABET[] =
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

std::string PathHasher::hashToBase62(uint64_t hash) {
    if (hash == 0) {
        return "0";
    }

    std::string result;
    result.reserve(11); // max length for 64-bit number in base62

    while (hash > 0) {
        result = BASE62_ALPHABET[hash % BASE62_BASE] + result;
        hash /= BASE62_BASE;
    }

    return result;
}

std::string PathHasher::hashPath(const std::string& path) {
    // Use std::hash for a simple, fast hash function
    std::hash<std::string> hasher;
    uint64_t hash = hasher(path);
    return hashToBase62(hash);
}
