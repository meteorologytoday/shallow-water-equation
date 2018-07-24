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

#include "configuration_shallow.hpp"
#include "fieldio.hpp"
#include "fftwfop.cpp" // template class must include its implementation
#include "vector_operation.cpp"




using namespace std;

#include "variables_declaration.hpp"

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

int main(int argc, char* args[]) {

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

	Q          = malloc_field_re();
	vort       = malloc_field_re();
	divg       = malloc_field_re();
	geop       = malloc_field_re();
	bg_vort    = malloc_field_re();

	u         = malloc_field_re();
	u2        = malloc_field_re();
	u_divg    = malloc_field_re();
	u_vort    = malloc_field_re();

    v         = malloc_field_re();
    v2        = malloc_field_re();
    v_divg    = malloc_field_re();
    v_vort    = malloc_field_re();


    K         = malloc_field_re();
    E         = malloc_field_re();

    absvort   = malloc_field_re();
    absvort_u = malloc_field_re();
    absvort_v = malloc_field_re();
    
    geop_u    = malloc_field_re();
    geop_v    = malloc_field_re();
    
    dvortdx   = malloc_field_re();
    dvortdy   = malloc_field_re();
    
    dvortdt   = malloc_field_re();
    ddivgdt   = malloc_field_re();
    dgeopdt   = malloc_field_re();
    
    workspace = malloc_field_re();


    // complex numbers
    
    vort_c0 = malloc_field_im();
    vort_c  = malloc_field_im(); 
    lvort_c = malloc_field_im();

    divg_c0 = malloc_field_im();
    divg_c  = malloc_field_im(); 
    ldivg_c = malloc_field_im();

    geop_c0 = malloc_field_im();
    geop_c  = malloc_field_im(); 
    lgeop_c = malloc_field_im();
    
    absvort_u_c  = malloc_field_im();
    absvort_v_c  = malloc_field_im();
    geop_u_c     = malloc_field_im();
    geop_v_c     = malloc_field_im();
    E_c          = malloc_field_im();
    chi_c        = malloc_field_im();
    psi_c        = malloc_field_im();
    tmp_c        = malloc_field_im();
    copy_for_c2r = malloc_field_im();

    dvortdt_c    = malloc_field_im();
    ddivgdt_c    = malloc_field_im();
    dgeopdt_c    = malloc_field_im();
    
    for(int k=0; k < 4; ++k) {
        rk4_vort_c[k] = malloc_field_im();
        rk4_divg_c[k] = malloc_field_im();
        rk4_geop_c[k] = malloc_field_im();
    }

    // plans
    p_fwd_vort       = crt_fwd_plan(Q,       Q_c);
    
    p_fwd_vort       = crt_fwd_plan(vort, vort_c);
    p_bwd_vort       = crt_bwd_plan(vort, vort_c);

    p_fwd_divg       = crt_fwd_plan(divg, divg_c);
    p_bwd_divg       = crt_bwd_plan(divg, divg_c);

    p_fwd_geop       = crt_fwd_plan(geop, geop_c);
    p_bwd_geop       = crt_bwd_plan(geop, geop_c);

    p_bwd_dvortdx    = crt_bwd_plan(dvortdx, tmp_c);
    p_bwd_dvortdy    = crt_bwd_plan(dvortdy, tmp_c);
    
    p_bwd_u          = crt_bwd_plan(u     , tmp_c);
    p_bwd_u_vort     = crt_bwd_plan(u_vort, tmp_c);
    p_bwd_u_divg     = crt_bwd_plan(u_divg, tmp_c);

    p_bwd_v          = crt_bwd_plan(v     , tmp_c);
    p_bwd_v_vort     = crt_bwd_plan(v_vort, tmp_c);
    p_bwd_v_divg     = crt_bwd_plan(v_divg, tmp_c);


    p_fwd_dvortdt    = crt_fwd_plan(dvortdt, dvortdt_c);
    p_fwd_ddivgdt    = crt_fwd_plan(ddivgdt, ddivgdt_c);
    p_fwd_dgeopdt    = crt_fwd_plan(dgeopdt, dgeopdt_c);
    
    p_fwd_absvort_u  = crt_fwd_plan(absvort_u, absvort_u_c);
    p_fwd_absvort_v  = crt_fwd_plan(absvort_v, absvort_v_c);

    p_fwd_geop_u     = crt_fwd_plan(geop_u, geop_u_c);
    p_fwd_geop_v     = crt_fwd_plan(geop_v, geop_v_c);

    p_fwd_E          = crt_fwd_plan(E, E_c);

		// read input
	Lx = LX;
	Ly = LY;
	dx = Lx / XPTS;
	dy = Ly / YPTS;

	sprintf(filename, "%s/%s", input.c_str(), init_file.c_str());
	readField(filename, vort, GRIDS);

	auto getDvortdt = [&](){
		
		// 1. Take laplacian (diffusion terms)
		fop.laplacian(vort_c, lvort_c); 
		fop.laplacian(divg_c, ldivg_c); 
		fop.laplacian(geop_c,       lgeop_c);

		// 2. Invert divergent and rotational flow
		fop.invertLaplacian(vort_c, psi_c);
		fop.invertLaplacian(divg_c, chi_c);
		
		// - rotation flow
		// -- calculate u_vort
		fop.grady(psi_c, tmp_c);
		fftwf_execute(p_bwd_u_vort); fftwf_backward_normalize(u_vort);
		for(int i=0; i<GRIDS;++i) { u_vort[i] = -u_vort[i]; }

		// -- calculate v_vort
		fop.gradx(psi_c, tmp_c);
		fftwf_execute(p_bwd_v_vort); fftwf_backward_normalize(v_vort);

		// - divergent flow
		// -- calculate u_div
		fop.gradx(chi_c, tmp_c);
		fftwf_execute(p_bwd_u_divg); fftwf_backward_normalize(u_divg);

		// -- calculate v_div
		fop.grady(chi_c, tmp_c);
		fftwf_execute(p_bwd_v_divg); fftwf_backward_normalize(v_divg);

		// - add together
		vop.add(u, u_vort, u_divg);
		vop.add(v, v_vort, v_divg);

		// 3. Calculate nonlinear multiplication in physical space

		// - invert vort, h
		fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
		fftwf_execute(p_bwd_geop); fftwf_backward_normalize(geop);

		// - vort to absvort
		vop.add(absvort, vort, bg_vort);

		// - calculate multiplication
		vop.mul(absvort_u, absvort, u);
		vop.mul(absvort_v, absvort, v);

		vop.mul(geop_u, geop,  u);
		vop.mul(geop_v, geop,  v);

		vop.mul(v2,   v,     v);
		vop.mul(u2,   u,     u);
		vop.add( K , u2,    v2);
		vop.add( E ,  K,  geop);
		
		// - forward transformation
		fftwf_execute(p_fwd_absvort_u);
		fftwf_execute(p_fwd_absvort_v);
		fftwf_execute(p_fwd_geop_u);
		fftwf_execute(p_fwd_geop_v);
		fftwf_execute(p_fwd_E);

		// - get spectral derivative

		// -- [vort]
		vop.set(dvortdt_c, 0.0f);

		// --- diffusion, Rayleigh friction
		vop.mul(tmp_c, lvort_c, NU);   vop.isub(dvortdt_c, tmp_c);
		vop.mul(tmp_c,  vort_c, MU);   vop.iadd(dvortdt_c, tmp_c);

		// --- rest
		fop.gradx(absvort_u_c, tmp_c); vop.isub(dvortdt_c, tmp_c);
		fop.grady(absvort_v_c, tmp_c); vop.isub(dvortdt_c, tmp_c);

		// -- [divg]
		vop.set(ddivgdt_c, 0.0f);

		// --- diffusion, Rayleigh friction
		vop.mul(tmp_c, ldivg_c, NU);   vop.isub(ddivgdt_c, tmp_c);
		vop.mul(tmp_c,  divg_c, MU);   vop.iadd(ddivgdt_c, tmp_c);

		// --- rest
		fop.gradx(absvort_v_c, tmp_c); vop.iadd(ddivgdt_c, tmp_c);
		fop.grady(absvort_u_c, tmp_c); vop.isub(ddivgdt_c, tmp_c);
		fop.laplacian(E_c, tmp_c);     vop.isub(ddivgdt_c, tmp_c);

		// -- [h]
		vop.set(dgeopdt_c, 0.0f);

		// --- Source
		// vop.isub(dgeop_c, Q_c);

		// --- diffusion
		vop.mul(tmp_c, lgeop_c, NU);   vop.isub(dgeopdt_c, tmp_c);

		// --- rest
		fop.gradx(geop_u_c, tmp_c);       vop.isub(dgeopdt_c, tmp_c);
		fop.grady(geop_v_c, tmp_c);       vop.isub(dgeopdt_c, tmp_c);

		/*
		#ifdef OUTPUT_GRAD_VORT
		if(debug) {
			sprintf(filename, "%s/dvortdx_step_%d.bin", output.c_str(), step);
			writeField(filename, dvortdx, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);
		}
		#endif
		*/

	};


	auto RK4_evolve_single = [&](fftwf_complex * updated, fftwf_complex * updated0, fftwf_complex * rk4_term, float coe) {
		for(int i=0; i<HALF_GRIDS; ++i) {
			updated[i][0] = updated0[i][0] + rk4_term[i][0] * coe;
			updated[i][1] = updated0[i][1] + rk4_term[i][1] * coe;
		}
	};

	auto RK4_last_step = [&](fftwf_complex * updated, fftwf_complex * updated0, fftwf_complex** rk4_term) {
		for(int i=0; i<HALF_GRIDS; ++i) {
			updated[i][0] = updated0[i][0] + (rk4_term[0][i][0] + 2.0f * rk4_term[1][i][0] + 2.0f * rk4_term[2][i][0] + rk4_term[3][i][0]) / 6.0f;
			updated[i][1] = updated0[i][1] + (rk4_term[0][i][1] + 2.0f * rk4_term[1][i][1] + 2.0f * rk4_term[2][i][1] + rk4_term[3][i][1]) / 6.0f;
		}
	};



	float RK4_step_coe[3] = {0.5f, 0.5f, 1.0f};

	auto RK4_run = [&]() {
		memcpy(vort_c0, vort_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup
		memcpy(divg_c0, divg_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup
		memcpy(geop_c0,       geop_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup

		for(int k = 0; k < 4 ; ++k) {
			getDvortdt();
			vop.mul(rk4_vort_c[k], dvortdt_c, dt); 
			vop.mul(rk4_divg_c[k], ddivgdt_c, dt) ;
			vop.mul(rk4_geop_c[k], dgeopdt_c, dt);
 
			if(k==3) { continue; }

			RK4_evolve_single(vort_c, vort_c0, rk4_vort_c[k], RK4_step_coe[k]);
			RK4_evolve_single(divg_c, divg_c0, rk4_divg_c[k], RK4_step_coe[k]);
			RK4_evolve_single(geop_c, geop_c0, rk4_geop_c[k], RK4_step_coe[k]);
		}

		RK4_last_step(vort_c, vort_c0, rk4_vort_c);
		RK4_last_step(divg_c, divg_c0, rk4_divg_c);
		RK4_last_step(geop_c, geop_c0, rk4_geop_c);

	};

	/*
	auto output_field_from_complex = [&](fftw_complex * data, char * filename) {
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, vort_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
			sprintf(filename, "%s/vort_step_%d.bin", output.c_str(), step);
			writeField(filename, vort, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(vort_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

		
	};
	*/

	printf("Program initialization complete.\n");

	// Preparation: 
	fftwf_execute(p_fwd_vort);
	fftwf_execute(p_fwd_divg);
	fftwf_execute(p_fwd_geop);
	fftwf_execute(p_fwd_Q);

	int record_flag = 0;
	for(int step = 0; step < total_steps; ++step) {

		printf("# Step %d, time = %.2f", step, step * dt);
		if( (record_flag = ((step % record_step) == 0)) ) { printf(", record now!");}
		printf("\n");

		if(record_flag) {

			// Output vort
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, vort_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
			sprintf(filename, "%s/vort_step_%d.bin", output.c_str(), step);
			writeField(filename, vort, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(vort_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

			
			// Output divg
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, divg_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_divg); fftwf_backward_normalize(divg);
			sprintf(filename, "%s/divg_step_%d.bin", output.c_str(), step);
			writeField(filename, divg, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(divg_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);
	
			// Output h
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, geop_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_geop); fftwf_backward_normalize(geop);
			sprintf(filename, "%s/geop_step_%d.bin", output.c_str(), step);
			writeField(filename, geop, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(geop_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

		}


		RK4_run();

		//fftwf_destroy_plan(p);
		//fftwf_free(in); fftwf_free(out);
	}

	fclose(log_fd);
	printf("Program ends. Congrats!\n");
	return 0;
}
