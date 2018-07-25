#include "fftwfop.hpp"

#ifndef FFTWF_OPERATION_CPP
#define FFTWF_OPERATION_CPP
template<int XPTS, int YPTS> fftwf_operation<XPTS,YPTS>::fftwf_operation(float Lx, float Ly) {
	this->gradx_coe = (float*) fftwf_malloc(sizeof(float) * XPTS);
	this->grady_coe = (float*) fftwf_malloc(sizeof(float) * this->HALF_YPTS);
	this->laplacian_coe_inverse = (float*) fftwf_malloc(sizeof(float) * this->HALF_GRIDS);
	this->laplacian_coe = (float*) fftwf_malloc(sizeof(float) * this->HALF_GRIDS);
	this->dealiasing_mask = (float*) fftwf_malloc(sizeof(float) * this->HALF_GRIDS);
	this->dealiase_xwavenumber = (int) ceil(((float)XPTS)/3.0);
	this->dealiase_ywavenumber = (int) ceil(((float)YPTS)/3.0);



    // For nonlinear multiplication
	this->dealiased_in1_c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * this->HALF_GRIDS);
	this->dealiased_in2_c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * this->HALF_GRIDS);
	this->dealiased_in1   = (float*) fftwf_malloc(sizeof(float) * this->GRIDS);
	this->dealiased_in2   = (float*) fftwf_malloc(sizeof(float) * this->GRIDS);

    this->p_bwd_in1 = fftwf_plan_dft_c2r_2d(XPTS, YPTS, this->dealiased_in1_c, this->dealiased_in1, FFTW_ESTIMATE);
    this->p_bwd_in2 = fftwf_plan_dft_c2r_2d(XPTS, YPTS, this->dealiased_in2_c, this->dealiased_in2, FFTW_ESTIMATE);
    this->p_fwd_in1 = fftwf_plan_dft_r2c_2d(XPTS, YPTS, this->dealiased_in1, this->dealiased_in1_c, FFTW_ESTIMATE);





	// Creating partial derivatives
	for(int i=0; i < this->HALF_XPTS; ++i) {
		this->gradx_coe[i] = TWOPI * ((float) i) / Lx;
	}
	for(int i=this->HALF_XPTS; i<XPTS; ++i) {
		this->gradx_coe[i] = - this->gradx_coe[this->reflectedXWavenumberIndex(i)];
	}

	for(int j=0; j < this->HALF_YPTS; ++j) {
		this->grady_coe[j] = TWOPI * ((float) j) / Ly;
	}

	#ifndef NDEBUG
	printf("gradx_coe list: \n");
	for(int i=0; i < XPTS; ++i) {
		printf("%f ", this->gradx_coe[i]);
	}
	printf("\ngrady_coe list:\n");
	for(int j=0; j < this->HALF_YPTS; ++j) {
		printf("%f ", grady_coe[j]);
	}
	printf("\n");
	#endif


	// Creating Laplacian coefficients
	for(int i=0; i<this->HALF_XPTS; ++i) {
		for(int j=0; j<this->HALF_YPTS; ++j) {
			this->laplacian_coe_inverse[this->HIDX(i,j)] = - (pow(this->gradx_coe[i],2) + pow(this->grady_coe[j],2));
			if(i == 0 && j == 0) {this->laplacian_coe_inverse[this->HIDX(i,j)] = 1.0f;} // this coe is special

			this->laplacian_coe[this->HIDX(i,j)] = - (pow(this->gradx_coe[i],2) + pow(this->grady_coe[j],2));
		}
	}

	for(int i=this->HALF_XPTS; i < XPTS; ++i) {
		for(int j=0; j<this->HALF_YPTS; ++j) {
			this->laplacian_coe_inverse[this->HIDX(i,j)] = this->laplacian_coe_inverse[this->R_HIDX(i,j)];
			this->laplacian_coe[this->HIDX(i,j)] = this->laplacian_coe[this->R_HIDX(i,j)];
		}
	}

	// Creating mask
	float generalized_wavenumber_square = pow(this->dealiase_xwavenumber,2) + pow(this->dealiase_ywavenumber,2);

	for(int i = 0; i < this->HALF_XPTS ; ++i) {
		for(int j = 0; j < this->HALF_YPTS; ++j) {
			this->dealiasing_mask[this->HIDX(i,j)] = (pow(i,2) + pow(j,2) >= generalized_wavenumber_square) ? 0.0f : 1.0f;
		}
	}
	for(int i = this->HALF_XPTS; i<XPTS ; ++i) {
		for(int j = 0; j < HALF_YPTS; ++j) {
			this->dealiasing_mask[this->HIDX(i,j)] = this->dealiasing_mask[this->R_HIDX(i,j)];
		}
	}

	#ifndef NDEBUG
	for(int i = 0; i < XPTS ; ++i) {
		for(int j = 0; j < HALF_YPTS; ++j) {
			printf("%d ", (int)(this->dealiasing_mask[this->HIDX(i,j)]));
		}
		printf("\n");
	}
	#endif

}

template<int XPTS, int YPTS> fftwf_operation<XPTS,YPTS>::~fftwf_operation() {
	fftwf_free(this->gradx_coe);
	fftwf_free(this->grady_coe);
	fftwf_free(this->laplacian_coe);
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::gradx(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<XPTS; ++i) {
		for(int j=0; j<this->HALF_YPTS; ++j) {
			out[this->HIDX(i,j)][0] = - in[this->HIDX(i,j)][1] * this->gradx_coe[i];
			out[this->HIDX(i,j)][1] =   in[this->HIDX(i,j)][0] * this->gradx_coe[i];
		}
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::grady(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<XPTS; ++i) {
		for(int j=0; j<this->HALF_YPTS; ++j) {
			out[this->HIDX(i,j)][0]   = - in[this->HIDX(i,j)][1] * this->grady_coe[j];
			out[this->HIDX(i,j)][1]   =   in[this->HIDX(i,j)][0] * this->grady_coe[j];
		}
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::laplacian(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<this->HALF_GRIDS; ++i) {
		out[i][0] = in[i][0] * this->laplacian_coe[i];
		out[i][1] = in[i][1] * this->laplacian_coe[i];
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::invertLaplacian(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<this->HALF_GRIDS; ++i) {
		out[i][0] = in[i][0] / this->laplacian_coe_inverse[i];
		out[i][1] = in[i][1] / this->laplacian_coe_inverse[i];
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::dealiase(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<this->HALF_GRIDS; ++i) {
		out[i][0] = in[i][0] * this->dealiasing_mask[i];
		out[i][1] = in[i][1] * this->dealiasing_mask[i];
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::pseudospectral_mul(fftwf_complex *in1, fftwf_complex *in2, fftwf_complex *out) {
    this->dealiase(in1, this->dealiased_in1_c);
    this->dealiase(in2, this->dealiased_in2_c);

    fftwf_execute(this->p_bwd_in1);
    fftwf_execute(this->p_bwd_in2);
  
   	for(int i=0; i<this->GRIDS; ++i) {
        //printf("%f, %f, %d, %d\n", this->dealiased_in1[i], this->dealiased_in2[i],  GRIDS2, this->GRIDS);
	    this->dealiased_in1[i] = this->dealiased_in1[i] * this->dealiased_in2[i] / this->GRIDS / this->GRIDS;

    }

    fftwf_execute(this->p_fwd_in1);
   
   	for(int i=0; i < HALF_GRIDS; ++i) {
		out[i][0] = this->dealiased_in1_c[i][0];
		out[i][1] = this->dealiased_in1_c[i][1];
        //printf("%f, %f\n", out[i][0], out[i][1]);

	}

}


#endif
