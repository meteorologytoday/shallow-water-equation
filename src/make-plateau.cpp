#include <cstdio>
#include <cstdlib>

#define _USE_MATH_DEFINES
#include <cmath>
#include <errno.h>

#include "configuration.hpp"
#include "fieldio.hpp"


int main() {
	
	float centery = LY / 2.0, plateau_wid = 50000.0, skirt_wid = 10000.0,
          geop_base = 9.8 * 10000.0, geop_diff = 9.8 * 1000.0;


	float* vort = (float*) malloc(sizeof(float) * GRIDS);
	float* divg = (float*) malloc(sizeof(float) * GRIDS);
	float* geop = (float*) malloc(sizeof(float) * GRIDS);

    float inner = plateau_wid / 2.0, outer = inner + skirt_wid;


	float x, y, s;
    printf("%f\n", s);
	for(int i=0; i<XPTS; ++i) {
		x = i * dx;
		for(int j=0; j<YPTS; ++j) {
			y = j * dy;

            s = (abs(y - centery) - inner) / skirt_wid;
            printf("%f\n", s);
			if(s <= 0) {
				geop[IDX(i,j)] = geop_base;
			} else if (s <= 1) {
				geop[IDX(i,j)] = geop_base + (2.0 * pow(s, 3) - 3.0 * pow(s, 2)) * geop_diff;
			} else {
				geop[IDX(i,j)] = geop_base - geop_diff;
			}
		}
	}

	std::string file;

    file = input + "/init_divg.bin";
	writeField(file.c_str(), divg, GRIDS);

	file = input + "/init_vort.bin";
	writeField(file.c_str(), vort, GRIDS);

	file = input + "/init_geop.bin";
	writeField(file.c_str(), geop, GRIDS);

	return 0;
}

