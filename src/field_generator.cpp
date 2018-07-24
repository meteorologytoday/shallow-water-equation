#define _USE_MATH_DEFINES
#include <cmath>

#include "configuration.hpp"

float dist(float x, float y, float cx, float cy) {
	return sqrtf(pow(x-cx,2.0) + pow(y-cy,2.0));
}

void addCakeKuo2004(float * data, float cx, float cy, float zeta_0, float scale_r) {

	float x, y, r;

	for(size_t j=0; j < XPTS; ++j) {

		y = j * DY;

		for(size_t i=0; i < YPTS; ++i) {

			x = i * DX;
			r = dist(x, y, cx, cy) / scale_r;

			if(r < 1) {
				data[IDX(i,j)] += zeta_0 * ( 1 - exp(- 30.0 / r * exp(1.0 / (r - 1.0))) );
			} 
		}
	}
}
