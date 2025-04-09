/**
 * File: database_writer.h
 * Created: 2022-12-09 14:53:33
 * Author: Milot Mirdita (milot@mirdita.de)
 */

#ifndef DATABASE_WRITER_H
#define DATABASE_WRITER_H
#include <cstdint>
#include <cstddef>


void* make_writer(const char *data_name, const char *index_name);
void free_writer(void *reader);

bool writer_append(void *reader, const char* data, size_t length, uint32_t key, const char* name);

#endif
