/**
 * File: database_reader.cpp
 * Created: 2023-02-10 17:04:07
 * Author: Milot Mirdita (milot@mirdita.de)
 */

#include "database_reader.h"
#include "utility.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <sys/stat.h>

struct reader_index_s {
    uint32_t id;
    int64_t length;
    int64_t offset;
};
typedef struct reader_index_s reader_index;

typedef std::vector<std::pair<std::string, uint32_t>> lookup_entry;

struct DBReader_s {
    reader_index* index;
    int64_t size;

    char* data;
    int64_t data_size;

    int dataMode;

    bool cache;

    lookup_entry* lookup;
};
typedef struct DBReader_s DBReader;

enum {
    SORT_BY_FIRST = 0,
    SORT_BY_SECOND = 1
};

ssize_t count_lines(char *data, ssize_t size);
struct compare_by_id {
    bool operator()(const reader_index &a, const reader_index &b) const {
        return a.id < b.id;
    }
};
bool read_index(DBReader *reader, char *data);
bool read_lookup(lookup_entry &lookup, char *data, ssize_t size, int sortMode);
DBReader* load_cache(const char *name);
bool save_cache(DBReader *reader, const char *name);

void* make_reader(const char *data_name, const char *index_name, int32_t data_mode) {
    char *data = NULL;
    ssize_t data_size = 0;
    if (data_mode & DB_READER_USE_DATA) {
        FILE* file = fopen(data_name, "r");
        if (file == NULL) {
            return NULL;
        }
        data = file_map(file, &data_size, 0);
        fclose(file);
    }

    // char cache_name[FILENAME_MAX];
    // if ((data_mode & DB_READER_NO_CACHE) == 0) {
    //     sprintf(cache_name, "%s.cache.%d", index_name, data_mode);

    //     struct stat st;
    //     if (stat(cache_name, &st) == 0) {
    //         DBReader* reader = load_cache(cache_name);
    //         reader->data = data;
    //         reader->data_size = data_size;
    //         reader->dataMode = data_mode;
    //         reader->cache = true;
    //         return (void*) reader;
    //     }
    // }


    FILE *file = fopen(index_name, "rb");
    if (file == NULL) {
        return NULL;
    }

    ssize_t index_size;
    char* index_data = file_map(file, &index_size, 0);
    DBReader* reader = (DBReader*)malloc(sizeof(DBReader));
    reader->size = count_lines(index_data, index_size);
	reader->index = (reader_index*)malloc(sizeof(reader_index) * reader->size);
	reader->data = data;
	reader->data_size = data_size;
	reader->dataMode = data_mode;
	reader->cache = false;
	if (!read_index(reader, index_data)) {
        free_reader(reader);
        return NULL;
    }
    file_unmap(index_data, (size_t)index_size);
    fclose(file);
    std::sort(reader->index, reader->index + reader->size, compare_by_id());

    reader->lookup = NULL;
    if (data_mode & (DB_READER_USE_LOOKUP) || (data_mode & DB_READER_USE_LOOKUP_REVERSE)) {
        std::string lookup_name(data_name);
        lookup_name = lookup_name + ".lookup";

        struct stat st;
        if (stat(lookup_name.c_str(), &st) == 0) {
            reader->lookup = new lookup_entry();
            reader->lookup->reserve(reader->size);
            FILE* file = fopen(lookup_name.c_str(), "rb");
            if (file == NULL) {
                free_reader(reader);
                return NULL;
            }
            ssize_t lookup_size;
            char *lookup_data = file_map(file, &lookup_size, 0);
            int sortMode = SORT_BY_FIRST;
            if (data_mode & DB_READER_USE_LOOKUP_REVERSE) {
                sortMode = SORT_BY_SECOND;
            }
            read_lookup(*(reader->lookup), lookup_data, lookup_size, sortMode);
            file_unmap(lookup_data, lookup_size);
            fclose(file);
        }
    }

    // if ((data_mode & DB_READER_NO_CACHE) == 0) {
    //     save_cache(reader, cache_name);
    // }

    return (void *)reader;
}

void free_reader(void *r) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL) {
        return;
    }

    if (reader->dataMode & DB_READER_USE_DATA) {
        file_unmap(reader->data, (size_t)(reader->data_size));
    }

    if (reader->cache) {
        file_unmap((char*)reader->index, (size_t)(reader->size) * sizeof(reader_index));
    } else {
        free(reader->index);
    }

    if (reader->lookup != NULL) {
        delete reader->lookup;
    }

    free(reader);
}

// ID is position in index and KEY pairs with the lookup name

int64_t reader_get_id(void *r, uint32_t key) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL) {
        return -1;
    }

    reader_index val;
    val.id = key;
    int64_t id = std::lower_bound(reader->index, reader->index + reader->size, val, compare_by_id()) - reader->index;
    if (id < reader->size && reader->index[id].id == key) {
        return id;
    } else {
        return -1;
    }
}

const char* reader_get_data(void *r, int64_t id) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL || id < 0 || id >= reader->size) {
        return NULL;
    }

    if (reader->index[id].offset >= reader->data_size) {
        return NULL;
    }

    return reader->data + reader->index[id].offset;
}

uint32_t reader_get_key(void *r, int64_t id) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL || id < 0 || id >= reader->size) {
        return -1;
    }
    return reader->index[id].id;
}

int64_t reader_get_length(void *r, int64_t id) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL || id < 0 || id >= reader->size) {
        return -1;
    }
    return reader->index[id].length;
}

int64_t reader_get_offset(void *r, int64_t id) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL || id < 0 || id >= reader->size) {
        return -1;
    }
    return reader->index[id].offset;
}

int64_t reader_get_size(void *r) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL) {
        return -1;
    }
    return reader->size;
}

ssize_t count_lines(char *data, ssize_t size) {
    size_t cnt = 0;
    for (ssize_t i = 0; i < size; ++i) {
        if (data[i] == '\n') {
            cnt++;
        }
    }
    return cnt;
}

size_t skipWhitespace(char * data) {
    size_t counter = 0;
    while ((data[counter] == ' ' || data[counter] == '\t') == true ) {
        counter++;
    }
    return counter;
}

size_t skipNoneWhitespace(char * data) {
    size_t counter = 0;
    while ((data[counter] == ' ' || data[counter] == '\t'
            || data[counter] == '\n' || data[counter] == '\0') == false ) {
        counter++;
    }
    return counter;
}

char* skipLine(char *data) {
     while (*data !='\n') {
        data++;
    }
    return (data+1);
}

size_t getWordsOfLine(char * data, char ** words, size_t maxElement ){
    size_t elementCounter = 0;
    while (*data != '\n' && *data != '\0'){
        data += skipWhitespace(data);
        words[elementCounter] = data;
        elementCounter++;
        if (elementCounter >= maxElement) {
            return elementCounter;
        }
        data += skipNoneWhitespace(data);
    }

    if(elementCounter < maxElement) {
        words[elementCounter] = data;
    }

    return elementCounter;
}

bool read_index(DBReader *reader, char *data) {
    bool status = true;
    int64_t i = 0;
    char *entry[255];
    while (i < reader->size) {
        const size_t columns = getWordsOfLine(data, entry, 255);

        if (columns > 3) {
            return false;
        }

        reader->index[i].id = (uint32_t)strtoul(entry[0], NULL, 10);
        int64_t offset = strtoull(entry[1], NULL, 10);
        int64_t length = strtoull(entry[2], NULL, 10);

        reader->index[i].length = length;

        if (reader->dataMode & DB_READER_USE_DATA) {
            reader->index[i].offset = offset;
        } else {
            reader->index[i].offset = 0;
        }

        i++;
        data = skipLine(data);
    }

    return status;
}

// compare_by_name
struct sort_by_first {
    bool operator()(const std::pair<std::string, uint32_t> &a, const std::pair<std::string, uint32_t> &b) const {
        return a.first.compare(b.first) <= 0;
    }
};
struct sort_by_second {
    bool operator()(const std::pair<std::string, uint32_t>& a, const std::pair<std::string, uint32_t>& b) const {
        return a.second <= b.second;
    }
};

struct compare_by_first {
    bool operator()(const std::pair<std::string, uint32_t> &lhs, const std::string &rhs) const {
        return  (lhs.first < rhs);
    }

    bool operator()(const std::string &lhs, const std::pair<std::string, uint32_t> &rhs) const {
        return  (lhs < rhs.first);
    }
};

struct compare_by_second {
    bool operator()(const std::pair<std::string, uint32_t>& lhs, const uint32_t& rhs) const {
        return  (lhs.second < rhs);
    }

    bool operator()(const uint32_t& lhs, const std::pair<std::string, uint32_t>& rhs) const {
        return  (lhs < rhs.second);
    }
};

bool read_lookup(lookup_entry &lookup, char *data, ssize_t size, int sortMode) {
    char *entry[3];
    ssize_t pos = 0;
    char* start = (char *) data;
    // size_t i = 0;
    while (pos < size) {
        const size_t columns = getWordsOfLine(data, entry, 3);
        if (columns < 3) {
            return false;
        }
        std::string name(entry[1], (entry[2] - entry[1]) - 1);
        uint32_t key = (uint32_t)strtoul(entry[0], NULL, 10);
        lookup.emplace_back(name, key);
        data = skipLine(data);
        pos = data - start;
        // i++;
    }
    if (sortMode == SORT_BY_FIRST) {
        std::stable_sort(lookup.begin(), lookup.end(), sort_by_first());
    } else {
        std::stable_sort(lookup.begin(), lookup.end(), sort_by_second());
    }
    return true;
}

uint32_t reader_lookup_entry(void* r, const char* name) {
    DBReader *reader = (DBReader*)r;
    if (reader == NULL || reader->lookup == NULL || reader->lookup->size() == 0) {
        return UINT32_MAX;
    }

    std::string name_str(name);
    lookup_entry::const_iterator it = std::lower_bound(reader->lookup->cbegin(), reader->lookup->cend(), name_str, compare_by_first());
    if (it != reader->lookup->cend() && it->first == name_str) {
        return it->second;
    }
    return UINT32_MAX;
}

const char* reader_lookup_name_alloc(void* r, uint32_t key) {
    DBReader* reader = (DBReader*)r;
    if (reader == NULL || reader->lookup == NULL || reader->lookup->size() == 0) {
        return "";
    }

    lookup_entry::const_iterator it = std::lower_bound(reader->lookup->cbegin(), reader->lookup->cend(), key, compare_by_second());
    if (it != reader->lookup->cend() && it->second == key) {
        return strdup(it->first.c_str());
    }
    return "";
}

DBReader* load_cache(const char *name) {
    FILE *file = fopen(name, "rb");
    if (file != NULL) {
        DBReader *reader = (DBReader*) malloc(sizeof(DBReader));
        ssize_t size;
        reader->index = (reader_index *) file_map(file, &size);
        reader->size = size / sizeof(reader_index);
        fclose(file);
        return reader;
    } else {
        return NULL;
    }
}

bool save_cache(DBReader *reader, const char *name) {
    FILE *file = fopen(name, "w+b");
    if (file != NULL) {
        fwrite(reader->index, sizeof(reader_index), (size_t)reader->size, file);
        fclose(file);
        return true;
    } else {
        return false;
    }
}
//
//int main(int argc, const char** argv) {
//    void* handle = make_reader("/Users/mirdita/tmp/pref", "/Users/mirdita/tmp/pref.index", 1);
//    int64_t id = reader_get_id(handle, 500);
//    printf("%lld\n", id);
//    printf("%s\n", reader_get_data(handle, id));
//    printf("%lld\n", reader_get_length(handle, id));
//    printf("%lld\n", reader_get_offset(handle, id));
//    free_reader(handle);
//}
