#include <fftw3.h>
#include <cmath>
#include <cassert>

#ifndef FFTWF_OPERATION_HPP
#define FFTWF_OPERATION_HPP
const float TWOPI = (acos(-1.0f) * 2.0f);

template<int XPTS, int YPTS> class fftwf_operation {
private:
	float *gradx_coe, *grady_coe, *laplacian_coe, *laplacian_coe_inverse, *dealiasing_mask;

    fftwf_complex *dealiased_in1_c, *dealiased_in2_c;
    float *dealiased_in1, *dealiased_in2;
    fftwf_plan p_bwd_in1, p_bwd_in2, p_fwd_in1;

	int dealiase_xwavenumber, dealiase_ywavenumber;
	int const HALF_XPTS = (int)(XPTS/2) + 1,
	    HALF_YPTS = (int)(YPTS/2) + 1,
		HALF_GRIDS = XPTS*HALF_YPTS,
		GRIDS = XPTS * YPTS,
        GRIDS2 = XPTS*XPTS*YPTS*YPTS;
public:
	fftwf_operation(float Lx, float Ly);
	~fftwf_operation();

	void gradx(fftwf_complex *in, fftwf_complex *out);
	void grady(fftwf_complex *in, fftwf_complex *out);
	void laplacian(fftwf_complex *in, fftwf_complex *out);
	void invertLaplacian(fftwf_complex *in, fftwf_complex *out);
	void dealiase(fftwf_complex *in, fftwf_complex *out);
    void pseudospectral_mul(fftwf_complex *in1, fftwf_complex *in2, fftwf_complex *out);

	inline int reflectedXWavenumberIndex(int i) {assert(i>=1 && "Input of ReflectedXWavenumberIndex must >= 1");return XPTS - i;};
	inline int HIDX  (int i, int j) {return HALF_YPTS*i + j;};
	inline int R_HIDX(int i, int j) {return HIDX(this->reflectedXWavenumberIndex(i),j);};
};


#endif
