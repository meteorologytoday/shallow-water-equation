/*
 * Author: Hsu, Tien-Yiao
 *
 * Description: 
 *
 * This program uses the psi output from main program to
 * find the extremum values of a binary file.
 *
 */

#include <cstdio>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "configuration.hpp"
#include "fieldio.hpp"

#define EXTREME_MIN 0
#define EXTREME_MAX 1


void trim(char * str) {
	char * p = str + strlen(str) - 1;
	while(p >= str) {
		if(strcmp(p,"\n") == 0) {
			*p = '\0';
		}
		--p;
	}
}

size_t find_max_pos(float * data, size_t sz) {
	size_t max_i = 0;
	for(size_t i = 1; i < sz; ++i) {
		if(data[i] > data[max_i]) {
			max_i = i;
		}
	}
	return max_i;
}

void find_min_n(float * data, size_t data_sz, float * result, size_t * result_pos, size_t result_sz) {

	size_t max_i;

	if(result_sz > data_sz) {
		fprintf(stderr, "Data size is %d, but you request %d numbers.\n", data_sz, result_sz);
		return;
	}

	for(size_t i=0; i < result_sz; ++i) {
		result[i] = data[i];
		result_pos[i] = i;
	}

	max_i = find_max_pos(result, result_sz);
	for(size_t i = result_sz; i < data_sz; ++i) {
		if(data[i] < result[max_i]) {
			result[max_i] = data[i];
			result_pos[max_i] = i;
			max_i = find_max_pos(result, result_sz);
		}
	}
}


int main(int argc, char* args[]) {
	fprintf(stderr, "Entering find_min program.\n");

	size_t min_n = 30;

	float * data = (float*) malloc(sizeof(float) * GRIDS);
	float * min = (float*) malloc(sizeof(float) * min_n);
	size_t * pos = (size_t*) malloc(sizeof(size_t) * min_n);

	// read input
	char filename[1024];
	size_t ptx, pty;
	while(fgets(filename, 1024, stdin) != NULL) {
		trim(filename);	
		readField(filename, data, GRIDS);
		fprintf(stderr, "File %s read.\n", filename);
		find_min_n(data, GRIDS, min, pos, min_n);

		for(size_t i = 0; i < min_n; ++i) {
			pty = pos[i] % YPTS;
			ptx = (size_t)(pos[i] / YPTS);
			fprintf(stdout, "%zu %zu %.5e\n", ptx, pty, min[i]);
		}

/*
		size_t written = fwrite(min, sizeof(float), min_n, stdout);
		if(written == min_n) {
				fprintf(stderr, "%d x %d bytes were written\n", written, sizeof(float));
		}
*/
	}

	fprintf(stderr, "find_min program ends. Congrats!\n");
	return 0;
}
