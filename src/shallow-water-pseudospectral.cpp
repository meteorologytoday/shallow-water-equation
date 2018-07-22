#define NDEBUG

#include <cstdio>

#define _USE_MATH_DEFINES
#include <cmath>

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <errno.h>
#include <assert.h>
#include <unistd.h> // getopt

#include "configuration.hpp"
#include "fieldio.hpp"

using namespace std;
using namespace VORT_SRC_READER;

float dx, dy, Lx, Ly;
float *vort, *u, *v, *dvortdx, *dvortdy, *dvortdt, *workspace, *vort_src;

string vort_src_filename = "";
int has_vort_src = 0;
enum RECIPE_TYPE recipe_type = EMPTY;

fftwf_plan p_fwd_vort,    p_bwd_vort,
		   p_bwd_dvortdx, p_bwd_dvortdy,
		   p_bwd_u,       p_bwd_v,
		   p_fwd_dvortdt, p_bwd_psi;

fftwf_complex *vort_c0, *vort_c, *lvort_c, *tmp_c, *psi_c, *rk1_c, *rk2_c, *rk3_c, *copy_for_c2r;

fftwf_operation<XPTS,YPTS> fop(LX, LY);

char filename[256];

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


int main(int argc, char* args[]) {
	char opt;

	while ((opt = getopt(argc, args, "I:O:i:s:f:")) != EOF) {
		switch(opt) {
			case 'I':
				input = optarg;
				break;
			case 'O':
				output = optarg;
				break;
			case 'i':
				init_file = optarg;
				break;
			case 's':
				vort_src_filename = optarg;
				recipe_type = SCRIPT;
				break;
			case 'f':
				vort_src_filename = optarg;
				recipe_type = FIFO;
				break;
		}
	}


	printf("##### Model setting #####\n");
	printf("Initial file          : %s \n", init_file.c_str());
	printf("Input folder          : %s \n", input.c_str());
	printf("Output folder         : %s \n", output.c_str());
	printf("Length X              : %.3f [m]\n", LX);
	printf("Length Y              : %.3f [m]\n", LY);
	printf("Spatial Resolution dx : %.3f [m]\n", dx);
	printf("Spatial Resolution dy : %.3f [m]\n", dy);
	printf("Time Resolution dt    : %.3f [s]\n", dt);
	printf("#########################\n\n\n");
	
	printf("Start project.\n");
	
	// open log
	FILE *log_fd = fopen("log", "w");
	if(log_fd == NULL) {
		perror("Open log file");
	}

	// initiate variables
	
	pv        = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	div       = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	h         = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	f         = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	u_rot     = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v_rot     = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	u_div     = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v_div     = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	
	u         = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v         = (float*) fftwf_malloc(sizeof(float) * GRIDS);


	// initializing plan
	p_fwd_vort       = fftwf_plan_dft_r2c_2d(XPTS, YPTS, vort, vort_c, FFTW_ESTIMATE);
	p_fwd_dvortdt    = fftwf_plan_dft_r2c_2d(XPTS, YPTS, dvortdt, tmp_c, FFTW_ESTIMATE);

	p_bwd_vort       = fftwf_plan_dft_c2r_2d(XPTS, YPTS, vort_c, vort, FFTW_ESTIMATE);
	p_bwd_dvortdx    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdx, FFTW_ESTIMATE);
	p_bwd_dvortdy    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdy, FFTW_ESTIMATE);
	p_bwd_u          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, u, FFTW_ESTIMATE);
	p_bwd_v          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, v, FFTW_ESTIMATE);

	p_bwd_psi        = fftwf_plan_dft_c2r_2d(XPTS, YPTS, psi_c, workspace, FFTW_ESTIMATE);


	dvortdx   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdy   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdt   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	workspace = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	vort_src  = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	// complex numbers
	vort_c0   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	vort_c    = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	lvort_c   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);  // laplacian vorticity complex
	tmp_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	psi_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk1_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk2_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk3_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	copy_for_c2r = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);

	// read input
	Lx = LX;
	Ly = LY;
	dx = Lx / XPTS;
	dy = Ly / YPTS;

	sprintf(filename, "%s/%s", input.c_str(), init_file.c_str());
	readField(filename, vort, GRIDS);

	auto getDvortdt = [&](bool debug, int step){


		
		// 1. Take laplacian (diffusion terms)
		fop.laplacian(vort_c, lvort_c);
		fop.laplacian(divg_c, ldivg_c);
		fop.laplacian(h_c,        h_c);

		// 2. Invert divergent and rotational flow
		fop.invertLaplacian(vort_c, psi_c);
		fop.invertLaplacian(divg_c, chi_c);
		
		
		// - calculate u_rot
		fop.grady(psi_c, tmp_c);
		fftwf_execute(p_bwd_u); fftwf_backward_normalize(u);
		for(int i=0; i<GRIDS;++i) { u[i] = -u[i]; }

		// - calculate v_rot
		fop.gradx(psi_c, tmp_c);
		fftwf_execute(p_bwd_v); fftwf_backward_normalize(v);


		// 3. 
		fftwf_execute(p_bwd_dvortdx); fftwf_backward_normalize(dvortdx);



		// step 03 take dvortdx, save as tmp_c
		fop.gradx(vort_c, tmp_c);

		// step 04
		fftwf_execute(p_bwd_dvortdx); fftwf_backward_normalize(dvortdx);

		#ifdef OUTPUT_GRAD_VORT
		if(debug) {
			sprintf(filename, "%s/dvortdx_step_%d.bin", output.c_str(), step);
			writeField(filename, dvortdx, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);
		}
		#endif

		// step 05 take dvortdy, save as tmp_c
		fop.grady(vort_c, tmp_c);

		// step 06
		fftwf_execute(p_bwd_dvortdy); fftwf_backward_normalize(dvortdy);

		#ifdef OUTPUT_GRAD_VORT
		if(debug) {
			sprintf(filename, "%s/dvortdy_step_%d.bin", output.c_str(), step);
			writeField(filename, dvortdy, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);
		}
		#endif

		// step 07 get psi_c
		fop.invertLaplacian(vort_c, psi_c);

		#ifdef OUTPUT_PSI

		if(debug) {
			 // backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, psi_c, sizeof(fftwf_complex) * HALF_GRIDS);
			fftwf_execute(p_bwd_psi); fftwf_backward_normalize(workspace);
			sprintf(filename, "%s/psi_step_%d.bin", output.c_str(), step);
			writeField(filename, workspace, GRIDS);
			fprintf(log_fd, "%s\n", filename);
			 // restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(psi_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);
		}
		
		#endif


		// step 08 get u_c
		fop.grady(psi_c, tmp_c);
		// step 09
		fftwf_execute(p_bwd_u); fftwf_backward_normalize(u);
		for(int i=0; i<GRIDS;++i) { u[i] = -u[i]; }

		#ifdef OUTPUT_WIND
		if(debug) {
			sprintf(filename, "%s/u_step_%d.bin", output.c_str(), step);
			writeField(filename, u, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);
		}
		#endif

		// step 10 get v_c
		fop.gradx(psi_c, tmp_c);
		// step 11
		fftwf_execute(p_bwd_v); fftwf_backward_normalize(v);

		#ifdef OUTPUT_WIND
		if(debug) {
			sprintf(filename, "%s/v_step_%d.bin", output.c_str(), step);
			writeField(filename, v, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);
		}
		#endif

		// step 12 get dvortdt
		for(int i=0; i<GRIDS;++i) {
			dvortdt[i] = - u[i] * dvortdx[i] - v[i] * dvortdy[i] + vort_src[i];
		}

		#ifdef OUTPUT_DVORTDT
		if(debug) {
			sprintf(filename, "%s/dvortdt_step_%d.bin", output.c_str(), step);
			writeField(filename, dvortdt, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);
		}
		#endif
		// step 13 get dvortdt_c and save in tmp_c
		fftwf_execute(p_fwd_dvortdt);

		// step ?? add laplacian term
		for(int i=0; i<HALF_GRIDS;++i) {
			tmp_c[i][0] += lvort_c[i][0] * NU;
			tmp_c[i][1] += lvort_c[i][1] * NU;
		}
	};

	auto evolve = [&](fftwf_complex *rk, float dt) {
		for(int i=0; i<HALF_GRIDS; ++i){
			vort_c[i][0] = vort_c0[i][0] + rk[i][0] * dt;
			vort_c[i][1] = vort_c0[i][1] + rk[i][1] * dt;
		}
	};

	printf("Initialization complete.\n");

	// step 01
	fftwf_execute(p_fwd_vort);

	// step 02 Everything is ready!!
	int record_flag = 0;
	for(int step = 0; step < total_steps; ++step) {

		printf("# Step %d, time = %.2f", step, step * dt);
		if( (record_flag = ((step % record_step) == 0)) ) { printf(", record now!");}
		printf("\n");

		if(record_flag) {

			sprintf(filename, "%s/vort_src_input_step_%d.bin", output.c_str(), step);
			writeField(filename, vort_src, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, vort_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
			sprintf(filename, "%s/vort_step_%d.bin", output.c_str(), step);
			writeField(filename, vort, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(vort_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);
		}



		// read source
		vs_reader.read(step * dt);

		memcpy(vort_c0, vort_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup

		for(int k = 0 ; k < 4; ++k) {

			getDvortdt(record_flag && (k==0), step);

			// step 14+15 dealiasing dvortdt_c and save to rk?_c)
			// DEPENDS ON RK?
			switch(k) {
				case 0:
					fop.dealiase(tmp_c, rk1_c);	evolve(rk1_c, dt / 2.0f);
					break;
				case 1:
					fop.dealiase(tmp_c, rk2_c);	evolve(rk2_c, dt / 2.0f);
					break;
				case 2:
					fop.dealiase(tmp_c, rk3_c);	evolve(rk3_c, dt);
					break;
				case 3:
					// Actually the variable rk4_c is not needed because this is the last call
					// we can simply replace rk4_c by tmp_c.
					fop.dealiase(tmp_c, tmp_c);

					// step 24: get new vort_c
					for(int i=0; i<HALF_GRIDS; ++i){
						vort_c[i][0] = vort_c0[i][0] + (rk1_c[i][0] + 2.0f * rk2_c[i][0] + 2.0f * rk3_c[i][0] + tmp_c[i][0]) * dt / 6.0f;
						vort_c[i][1] = vort_c0[i][1] + (rk1_c[i][1] + 2.0f * rk2_c[i][1] + 2.0f * rk3_c[i][1] + tmp_c[i][1]) * dt / 6.0f;
					}

					break;

			}
		}



		//fftwf_destroy_plan(p);
		//fftwf_free(in); fftwf_free(out);
	}

	fclose(log_fd);
	printf("Program ends. Congrats!\n");
	return 0;
}
