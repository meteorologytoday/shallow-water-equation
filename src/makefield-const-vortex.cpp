#include <cstdio>
#include <cstdlib>

#define _USE_MATH_DEFINES
#include <cmath>
#include <errno.h>

#include "configuration.hpp"
#include "fieldio.hpp"


int main() {
	
	float centerx = LX / 2.0, centery = LY / 2.0, r_bound = 6000.0, zeta0 = 2e-5;
	float dx = LX / XPTS, dy = LY / YPTS;
	float r;

	float* vort = (float*) malloc(sizeof(float) * GRIDS);

	auto radius = [centerx, centery](float x, float y) -> float {
		return sqrtf(pow(x-centerx,2) + pow(y-centery,2));
	};

	float x, y;
	for(int i=0; i<XPTS; ++i) {
		x = i * dx;
		for(int j=0; j<YPTS; ++j) {
			y = j * dy;
			r = radius(x,y);

			if(r <= r_bound) {
				vort[IDX(i,j)] = zeta0;
			} else {
				vort[IDX(i,j)] = 0;
			}
		}
	}

	std::string file = input + "/" + init_file;
	writeField(file.c_str(), vort, GRIDS);

	return 0;
}

