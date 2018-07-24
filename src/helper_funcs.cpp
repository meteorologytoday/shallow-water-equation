#include <fftw3.h>
void fftwf_backward_normalize(float *data) {
	for(int i=0; i < GRIDS; ++i) {
		data[i] /= GRIDS;
	}
}

float sumSqr(fftwf_complex *c) {
	float strength = 0;
	for(int i=0; i < HALF_GRIDS;++i){
		strength += pow(c[i][0],2) + pow(c[i][1],2);
	}
	return strength;
}

void print_spectrum(fftwf_complex *in) {
	for(int i=0; i < XPTS;++i){
		for(int j=0; j < HALF_YPTS;++j){
			printf("(%+6.2f, %+6.2f) ", in[HIDX(i,j)][0], in[HIDX(i,j)][1]);
		}
		printf("\n");
	}
}

void print_error(char * str) {
	printf("Error: %s\n", str);
}

float * malloc_field_re() {
	return (float*) fftwf_malloc(sizeof(float) * GRIDS);
}


fftwf_complex * malloc_field_im() {
    return (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
}

fftwf_plan crt_fwd_plan(float * re, fftwf_complex * im) {
    return fftwf_plan_dft_r2c_2d(XPTS, YPTS, re, im, FFTW_ESTIMATE);
}

fftwf_plan crt_bwd_plan(float * re, fftwf_complex * im) {
    return fftwf_plan_dft_c2r_2d(XPTS, YPTS, im, re, FFTW_ESTIMATE);
}


