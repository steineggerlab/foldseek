/**
 * File: database_writer.cpp
 * Created: 2022-12-09 14:53:34
 * Author: Milot Mirdita (milot@mirdita.de)
 */

#include "database_writer.h"
#include "database_reader.h"

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>

struct writer_index_s {
    uint32_t id;
    int64_t length;
    int64_t offset;
    uint32_t name_index;
};

typedef struct writer_index_s writer_index;

struct DatabaseWriter {
    FILE* data;
    FILE* index;
    FILE* lookup;
    writer_index* entries;
    std::vector<std::string> names;
    uint64_t size;
    uint64_t capacity;
    bool is_sorted;
};

void* make_writer(const char *data_name, const char *index_name) {
    DatabaseWriter* writer = new DatabaseWriter;
    writer->data = fopen(data_name, "wb");
    writer->index = fopen(index_name, "w");
    std::string lookup_name = std::string(data_name) + ".lookup";
    writer->lookup = fopen(lookup_name.c_str(), "w");
    writer->entries = (writer_index*)malloc(1000 * sizeof(writer_index));
    writer->size = 0;
    writer->capacity = 1000;
    writer->is_sorted = 1;
    // std::string source_name = std::string(data_name) + ".source";
    // FILE* source = fopen(source_name.c_str(), "w");
    // fprintf(source, "0\t%s", data_name);
    // fclose(source);
    std::string dbtype_name = std::string(data_name) + ".dbtype";
    FILE* dbtype = fopen(dbtype_name.c_str(), "w");
    // generic dbtype
    int type = 12;
    fwrite(&type, sizeof(int), 1, dbtype);
    fclose(dbtype);
    return writer;
}

void free_writer(void *writer) {
    DatabaseWriter* w = (DatabaseWriter*)writer;
    if (w->is_sorted == false) {
        std::stable_sort(w->entries, w->entries + w->size, [](const writer_index& a, const writer_index& b) { return a.id < b.id; });
    }
    for (uint64_t i = 0; i < w->size; ++i) {
        fprintf(w->index, "%d\t%llu\t%d\n", w->entries[i].id, w->entries[i].offset, (uint32_t)w->entries[i].length);
        fprintf(w->lookup, "%d\t%s\t0\n", w->entries[i].id, w->names[w->entries[i].name_index].c_str());
    }
    fclose(w->index);
    fclose(w->lookup);
    free(w->entries);
    fclose(w->data);
    delete w;
}

bool writer_append(void *writer, const char* data, size_t length, uint32_t key, const char* name) {
    DatabaseWriter* w = (DatabaseWriter*)writer;
    int64_t offset = ftell(w->data);
    size_t res = fwrite(data, 1, length, w->data);
    if (res != length) {
        return false;
    }
    writer_index entry;
    entry.id = key;
    entry.length = length;
    entry.offset = offset;
    w->names.push_back(name);
    entry.name_index = w->names.size() - 1;
    if (w->size == w->capacity) {
        w->capacity *= 2;
        w->entries = (writer_index*)realloc(w->entries, w->capacity * sizeof(writer_index));
    }
    w->entries[w->size] = entry;
    w->is_sorted = w->is_sorted && (w->size <= 1 || w->entries[w->size - 1].id < key);
    w->size++;
    return true;
}