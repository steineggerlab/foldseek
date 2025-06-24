/**
 * File: database_reader.h
 * Created: 2022-12-08 00:18:37
 * Author: Milot Mirdita (milot@mirdita.de)
 */

#ifndef DATABASE_READER_H
#define DATABASE_READER_H
#include <cstdint>

static const int DB_READER_USE_DATA   = 1u << 0;
static const int DB_READER_NO_CACHE   = 1u << 1;
static const int DB_READER_USE_LOOKUP = 1u << 2;
static const int DB_READER_USE_LOOKUP_REVERSE = 1u << 3;


void* make_reader(const char *data_name, const char *index_name, int32_t data_mode);
void free_reader(void *reader);

int64_t reader_get_id(void *reader, uint32_t key);
const char* reader_get_data(void *reader, int64_t id);
uint32_t reader_get_key(void *reader, int64_t id);
int64_t reader_get_length(void *reader, int64_t id);
int64_t reader_get_offset(void *reader, int64_t id);
int64_t reader_get_size(void *r);
uint32_t reader_lookup_entry(void* r, const char* name);
const char* reader_lookup_name_alloc(void* r, uint32_t key);

#endif
