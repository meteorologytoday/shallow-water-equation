#include <cstdio>
#include <cstdlib>

#define _USE_MATH_DEFINES
#include <cmath>
#include <errno.h>

#include "configuration.hpp"
#include "fieldio.hpp"


int main() {
	
	float centerx = LX / 2.0, centery = LY / 2.0, epsilon = 0.7, lambda = 2.0, zeta0 = .005f, r_i = 30000.0, r_o = 60000.0;
	float dx = LX / XPTS, dy = LY / YPTS;
	float r, r_i_alpha, r_o_alpha, r_prime;

	float* vort = (float*) malloc(sizeof(float) * GRIDS);

	auto radius = [centerx, centery](float x, float y) -> float {
		return sqrtf(pow(x-centerx,2) + pow(y-centery,2));
	};
	auto alpha = [centerx, centery, epsilon, radius](float x, float y) -> float {
		float c;
		float r = radius(x,y);
		if(r == 0.0f) {
			c = 0;
		} else {
			c = (y-centery) / radius(x,y);
		}
		return sqrtf((1.0 - pow(epsilon,2)) / (1.0 - pow(epsilon*( sin( asin(c) + 3.1415/4.0 ) ),2)));
	};

	float x, y;
	for(int i=0; i<XPTS; ++i) {
		x = i * dx;
		for(int j=0; j<YPTS; ++j) {
			y = j * dy;
			r = radius(x,y);
			r_i_alpha = r_i * alpha(x,y);
			r_o_alpha = r_o * alpha(x,y);

			if(r <= r_i_alpha) {
				vort[IDX(i,j)] = zeta0;
			} else if (r <= r_o_alpha) {
				r_prime = (r - r_i_alpha)/(r_o_alpha - r_i_alpha);
				vort[IDX(i,j)] = zeta0 * (1.0 - exp( - lambda / r_prime * exp(1.0 / (r_prime - 1))));
			} else {
				vort[IDX(i,j)] = 0;
			}
		}
	}

	std::string file = input + "/" + init_file;
	writeField(file.c_str(), vort, GRIDS);

	return 0;
}

