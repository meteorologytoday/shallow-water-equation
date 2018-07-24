#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>

#include "configuration.hpp"
#include "fieldio.hpp"
#include "field_generator.cpp"

// Two vortices (binary)
// 
// \Delta = 10 km
//
// vortex 1 (small, intense):
// 
//   >> R_1 = 10 km
//   >> zeta_1 = 1.5e-2 s^-1
//
// vortex 2 (large, weak):
//
//   >> R_2 = 30 km
//   >> zeta_2 = 3e-3 s^-1
//
// So, this imples center separation = 10 km + 10 km + 30 km = 50 km
//
//
//


int main() {
	
	float * vort = (float *) malloc(GRIDS * sizeof(float));
	

	float duration = 3600.0 * 3.0;
	float beg_time = 3600.0 * 2.0;
	float end_time = beg_time + duration;

	size_t beg_step = (size_t) (beg_time / dt);
	size_t end_step = (size_t) (end_time / dt);

	char flag;
	for(size_t step = 1 ; step < total_steps; ++step) {
		if(step == beg_step) {

			addCakeKuo2004(vort, LX/2.0 + 50000.0, LY/2.0, 3e-3 / duration, 30000.0);

			flag = (char) 1;
			fwrite(&flag, sizeof(char), 1, stdout);
			fwrite(vort, sizeof(float), GRIDS, stdout);
		
		} else if(step == end_step) {
			memset(vort, 0, GRIDS * sizeof(float) );
			flag = (char) 1;
			fwrite(&flag, sizeof(char), 1, stdout);
			fwrite(vort, sizeof(float), GRIDS, stdout);
		} else { 
			flag = (char) 0;
			fwrite(&flag, sizeof(char), 1, stdout);
		}
	}
	
	fprintf(stderr, "###### input program ends ######\n");

	return 0;
}

