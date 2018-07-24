#include <cstdio>
#include <cstdlib>

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

	// vortex 1
	addCakeKuo2004(vort, LX/2.0, LY/2.0, 1.5e-2, 10000.0);

	// vortex 2
//	addCakeKuo2004(vort, LX/2.0 + 50000.0, LY/2.0, 3e-3, 30000.0);
	
	std::string file = input + "/" + init_file;
	writeField(file.c_str(), vort, GRIDS);

	return 0;
}

