#ifndef PROST_T5_H
#define PROST_T5_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ProstT5 ProstT5;

void prostt5_free(ProstT5* ptr);
ProstT5* prostt5_load(const char* base_path, bool profile, bool cpu, bool cache);
const char* prostt5_predict(ProstT5* ptr, const char* sequence);
const char* prostt5_predict_slice(ProstT5* ptr, const char* sequence, size_t length);
void prostt5_free_cstring(const char* str);

#ifdef __cplusplus
}
#endif

#endif // PROST_T5_H

