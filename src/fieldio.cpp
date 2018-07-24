#include "fieldio.hpp"
#include <cstdio>
#include <cerrno>

#ifndef FIELDIO_CPP
#define FIELDIO_CPP
void writeField(const char * filename, float *data, size_t len) {

	FILE * file = fopen(filename, "wb");
	int flg;
	errno = 0;
	if((flg = fwrite(data, sizeof(float), len, file)) < 0) {
		perror("Write field.");
		//printf("ERROR! %d\n", flg);
	}
	fclose(file);

	fprintf(stderr, "Output %s\n", filename);
}

void readField(const char * filename, float *data, size_t len) {

	FILE * file = fopen(filename, "rb");
	int flg;
	errno = 0;
	if((flg = fread(data, sizeof(float), len, file)) < 0) {
		perror("Read field.");
		// printf("ERROR! %d\n", flg);
	}
	fclose(file);

	fprintf(stderr, "%d bytes read: %s\n", flg, filename);
}

#endif
