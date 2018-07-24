#include <cstddef>

#ifndef FIELDIO_H
#define FIELDIO_H
void writeField(const char * filename, float *data, size_t len);
void readField(const char * filename, float *data, size_t len);
#endif
