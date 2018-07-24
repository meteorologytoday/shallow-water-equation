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
#include "fftwfop.cpp" // template class must include its implementation
#include "vector_operation.cpp"

using namespace std;

#include "variables_declaration.cpp"
#include "helper_funcs.cpp"

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
	printf("GRIDS                 : %d\n", GRIDS);
	printf("HALF_GRIDS            : %d\n", HALF_GRIDS);
	printf("#########################\n\n\n");
	
	printf("Start project.\n");
	
	// open log
	FILE *log_fd = fopen("log", "w");
	if(log_fd == NULL) {
		perror("Open log file");
	}

	// initiate variables

    #include "variables_initialization.cpp"

    printf("Reading input... ");

	sprintf(filename, "%s/init_vort.bin", input.c_str());
	readField(filename, vort, GRIDS);

	//sprintf(filename, "%s/init_divg.bin", input.c_str());
	//readField(filename, vort, GRIDS);

	//sprintf(filename, "%s/init_geop.bin", input.c_str());
	//readField(filename, vort, GRIDS);
    printf("done.\n");


	auto getDvortdt = [&](){
	    //printf("Doing linear terms..."); fflush(stdout);	
		// 1. Take laplacian (diffusion terms)
		fop.laplacian(vort_c, lvort_c); 
		//fop.laplacian(divg_c, ldivg_c); 
		//fop.laplacian(geop_c, lgeop_c);

		// 2. Invert divergent and rotational flow
		fop.invertLaplacian(vort_c, psi_c);
		//fop.invertLaplacian(divg_c, chi_c);
		
	    //printf("done1.\n"); fflush(stdout);	
		// - rotation flow
		// -- calculate u_vort
		fop.grady(psi_c, tmp_c);
	    //printf("done2.\n"); fflush(stdout);	
        cvop.mul(u_c, tmp_c, -1.0);
	    //printf("done3.\n"); fflush(stdout);	
		// -- calculate v_vort
		fop.gradx(psi_c, v_c);

		// - divergent flow
		// -- calculate u_div
		//fop.gradx(chi_c, tmp_c);
        //cvop.iadd(u_c, tmp_c);

		// -- calculate v_div
		//fop.grady(chi_c, tmp_c);
        //cvop.iadd(v_c, tmp_c);

	    //printf("done.\n"); fflush(stdout);	

	    //printf("Doing nonlinear terms..."); fflush(stdout);	
		// 3. Calculate nonlinear multiplication in physical space

        cvop.add(absvort_c, vort_c, bg_vort_c);
	    //printf("done3.1.\n"); fflush(stdout);	
        fop.pseudospectral_mul(absvort_c, u_c, absvort_u_c);
	    //printf("done3.2.\n"); fflush(stdout);	
        fop.pseudospectral_mul(absvort_c, v_c, absvort_v_c);

	    //printf("done4.\n"); fflush(stdout);	
        //fop.pseudospectral_mul(geop_c, u_c, geop_u_c);
        //fop.pseudospectral_mul(geop_c, v_c, geop_v_c);
        
	    //printf("done5.\n"); fflush(stdout);	
        //fop.pseudospectral_mul(u_c, u_c, u2_c);
	    //printf("done5.5\n"); fflush(stdout);	
        //fop.pseudospectral_mul(v_c, v_c, v2_c);
	    //printf("done6.\n"); fflush(stdout);	
        //cvop.add(E_c, u2_c, v2_c);
        //cvop.iadd(E_c, geop_c);
	    //printf("done7.\n"); fflush(stdout);	
		
        // -- [vort]
		cvop.set(dvortdt_c, 0.0f);

		// --- diffusion, Rayleigh friction
		cvop.mul(tmp_c, lvort_c, NU);   cvop.isub(dvortdt_c, tmp_c);
		cvop.mul(tmp_c,  vort_c, MU);   cvop.iadd(dvortdt_c, tmp_c);

		// --- rest
		fop.gradx(absvort_u_c, tmp_c); cvop.isub(dvortdt_c, tmp_c);
		fop.grady(absvort_v_c, tmp_c); cvop.isub(dvortdt_c, tmp_c);

		// -- [divg]
		//cvop.set(ddivgdt_c, 0.0f);

		// --- diffusion, Rayleigh friction
		//cvop.mul(tmp_c, ldivg_c, NU);   cvop.isub(ddivgdt_c, tmp_c);
		//cvop.mul(tmp_c,  divg_c, MU);   cvop.iadd(ddivgdt_c, tmp_c);

		// --- rest
		//fop.gradx(absvort_v_c, tmp_c); cvop.iadd(ddivgdt_c, tmp_c);
		//fop.grady(absvort_u_c, tmp_c); cvop.isub(ddivgdt_c, tmp_c);
		//fop.laplacian(E_c, tmp_c);     cvop.isub(ddivgdt_c, tmp_c);

		// -- [geop]
		cvop.set(dgeopdt_c, 0.0f);

		// --- Source
		// cvop.isub(dgeop_c, Q_c);

		// --- diffusion
		//cvop.mul(tmp_c, lgeop_c, NU);   cvop.isub(dgeopdt_c, tmp_c);

		// --- rest
		//fop.gradx(geop_u_c, tmp_c);       cvop.isub(dgeopdt_c, tmp_c);
		//fop.grady(geop_v_c, tmp_c);       cvop.isub(dgeopdt_c, tmp_c);

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
	    printf("RK4 run\n"); fflush(stdout);
		memcpy(vort_c0, vort_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup
		memcpy(divg_c0, divg_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup
		memcpy(geop_c0,       geop_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup

		for(int k = 0; k < 4 ; ++k) {
			getDvortdt();
			cvop.mul(rk4_vort_c[k], dvortdt_c, dt); 
			cvop.mul(rk4_divg_c[k], ddivgdt_c, dt) ;
			cvop.mul(rk4_geop_c[k], dgeopdt_c, dt);
 
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
	printf("[1] Fwd Q ... "); fflush(stdout);
	fftwf_execute(p_fwd_Q);
	printf("done\n");

	int record_flag = 0;
	for(int step = 0; step < total_steps; ++step) {

		printf("# Step %d, time = %.2f", step, step * dt);
		if( (record_flag = ((step % record_step) == 0)) ) { printf(", record now!");}
		printf("\n");

        // NOW FOR TESTING, SET DIVERGENT = 0
        cvop.set(divg_c, 0.0f);


		if(record_flag) {

			// Output vort
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, vort_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
			sprintf(filename, "%s/vort_step_%04d.bin", output.c_str(), step);
			writeField(filename, vort, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(vort_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

			
			// Output divg
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, divg_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_divg); fftwf_backward_normalize(divg);
			sprintf(filename, "%s/divg_step_%04d.bin", output.c_str(), step);
			writeField(filename, divg, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(divg_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);
	
			// Output h
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, geop_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_geop); fftwf_backward_normalize(geop);
			sprintf(filename, "%s/geop_step_%04d.bin", output.c_str(), step);
			writeField(filename, geop, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(geop_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

		}


		RK4_run();

	}

	fclose(log_fd);
	printf("Program ends. Congrats!\n");
	return 0;
}
