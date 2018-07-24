#define DEBUG_MODE
#define X_INCLUDE_IO

#include "xdmlab"
#include <iostream>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <ctime>


#ifdef DEBUG_MODE
	#define __DBG__ 
#else
	#define __DBG__ //
#endif


using namespace X;
using namespace std;


typedef float T;
const xsize xpts = 512;
const xsize ypts = 512;
const T Lx = 300000.0f;
const T Ly = 300000.0f;
const T nu_vort = 0.0;//0.1;
const T g0 = 10;

T dx, dy;
const char * datafolder = "./output";


typedef DataIndex2<T,xpts,ypts> ModelData;
typedef Enumerator2<xpts, ypts> EnuData;

typedef DataIndex1<T,xpts> GridX;
typedef DataIndex1<T,ypts> GridY;


GridX x;
GridY y;
ModelData psi, vort, vort_tmp, vort_tmp2, k1_vort, k2_vort, k3_vort, k4_vort, u, v, tmp, gradx_vort, grady_vort, lap_vort, lap_theta, bc_filter;

const xsize step = 3600;
const xsize record_step = 10;
const float dt = 3.0f;
const xsize max_iter = 10000;
const T relative_err = 1e-4;

void init();
void gradx(ModelData&, ModelData&, T);
void grady(ModelData&, ModelData&, T);
int jacobi_relaxation(ModelData&, ModelData&, ModelData&, T, T, T, xsize);
void lap(ModelData&, ModelData&, T, T);
int output_data(ModelData&, const char *);


void init() {
	X::Util::dLinspace1(x,0.0f,Lx);
	X::Util::dLinspace1(y,0.0f,Ly);
	psi = 0.0f;
	vort = 0.0f;
	dx = x(1) - x(0);
	dy = y(1) - y(0);

	auto add_vort = [&] (ModelData& target_vort, T zeta0, T r_in, T r_out, T centx, T centy) {
		T r, tmp;
		EnuData::each_index([&](xsize i, xsize j){
			r = sqrt(pow(x(i) - centx,2) + pow(y(j) - centy, 2));
			if(r < r_in) {
				target_vort(i,j) += zeta0;
			} else if(r <= r_out) {
				tmp = (r - r_in)/(r_out - r_in);
				target_vort(i,j) += zeta0 * (2*tmp*tmp*tmp - 3*tmp*tmp + 1);
			} else {
				target_vort(i,j) += 0;
			}
		});
	};

	// create a sin(2 pi x/L) * sin(4 pi y / L) field
	T centx1 = 150000.0f, centy1 = 150000.0f;
	T centx2 = 195000.0f, centy2 = 150000.0f;
	T r_in_1 = 10000.f, r_out_1 = 12000.f;
	T r_in_2 = 25000.f, r_out_2 = 29000.f;
	T zeta0_1 = 0.02;
	T zeta0_2 = 0.003;

	add_vort(vort, zeta0_1, r_in_1, r_out_1, centx1, centy1);
	add_vort(vort, zeta0_2, r_in_2, r_out_2, centx2, centy2);

	EnuData::each_index([&](xsize i, xsize j){
		bc_filter(i,j) = ( i == 0 || j ==0 || i == xpts-1 || j == ypts-1 ) ? 0.0f : 1.0f;
	});


}

void gradx(ModelData& input, ModelData& output, T dx) {
	T d2x = 2*dx;
	for(xsize i=0; i<xpts; ++i) {
		for(xsize j=0; j<ypts; ++j) {
			if(i == 0) {
				output(i,j) = (input(i+1,j) - input(i,j)) / dx;
			} else if(i == xpts-1) {
				output(i,j) = (input(i,j) - input(i-1,j)) / dx;
			} else {
				output(i,j) = (input(i+1, j) - input(i-1,j)) / d2x;
			}
		}
	}
}

void grady(ModelData& input, ModelData& output, T dy) {
	T d2y = 2*dy;
	for(xsize i=0; i<xpts; ++i) {
		for(xsize j=0; j<ypts; ++j) {
			if(j == 0) {
				output(i,j) = (input(i,j+1) - input(i,j)) / dy;
			} else if(j == ypts-1) {
				output(i,j) = (input(i,j) - input(i,j-1)) / dy;
			} else {
				output(i,j) = (input(i, j+1) - input(i,j-1)) / d2y;
			}
		}
	}
}

void print_add(ModelData& d) {
	cout << "In Index: " << d.getTargetDataPointer() << ", in Data:" << d.p() << endl; 
}

int jacobi_relaxation(ModelData& input, ModelData& output, ModelData& tmp, T dx, T dy, T strategy, T strategy_value, xsize max_iter) {
	/*
	 * Boundary condition is specify on the margin of input
	 */
	T err = 0, err_before, r;
	T gamma = 1 / (2 * ( 1/(dx*dx) + 1/(dy*dy) ) );
	for(xsize iter = 0; iter < max_iter; ++iter) {
		err = 0;
		// calculate nth step zeta field
		lap(output, tmp, dx, dy);

		// tune
		for(xsize i = 1; i < xpts-1; ++i) {
			for(xsize j = 1; j < ypts-1; ++j) {
				r = tmp(i,j) - input(i,j);
				err += r * r;
				output(i,j) += gamma * r;
			}
		}
		err = sqrt(err/ ((xpts-2)*(ypts-2)));
		if(strategy == 0 && err < strategy_value) {
			return 0;
		} else if(strategy == 1) {
			if(iter > 0 && ((err_before - err) / err < strategy_value)) {
				cout << "iter: " << iter << "; err:" << err << endl;
				return 0;
			}

			err_before = err;
		}
	}
	
	return -1;
}

void lap(ModelData& input, ModelData& output, T dx, T dy) {
	/*
	 * Points must be larger than 3 in each direction.
	 */

	T dxdx = dx * dx;
	T dydy = dy * dy;
	xsize centx, centy;
	for(xsize i = 0; i < xpts; ++i) {
		for(xsize j = 0; j < ypts; ++j) {
			centx = i; centy = j;
			if(i == 0)
				++centx;
			else if(i == xpts-1)
				--centx;

			if(j == 0)
				++centy;
			else if(j == ypts-1)
				--centy;

			output(i,j) = (input(centx+1,centy) + input(centx-1,centy) - 2*input(centx,centy)) / dxdx \
						  +  (input(centx,centy+1) + input(centx,centy-1) - 2*input(centx,centy)) / dydy;
		}
	}

}

int output_data(ModelData& dat, const char * filename) {
	FILE * file = fopen(filename,"wb");
	T *buf = new T[xpts];
	for(xsize j = 0; j < ypts; ++j) {
		for(xsize i=0; i < xpts; ++i) {
			buf[i] = dat(i,j);
		}
		fwrite(buf, sizeof(T), xpts, file);
	}

	fclose(file);
	return 0;
}

int main() {

	init();

	auto ArakawaJacobian = [] (ModelData& p, ModelData& q, ModelData& JJ) {

		T p1,p2,p3,p4,p5,p6,p7,p8,\
		  q1,q2,q3,q4,q5,q6,q7,q8;
		int sh_le, sh_ri, sh_up, sh_do;
		EnuData::each_index( [&] (xsize i, xsize j) {
			sh_le = (i == 0     ) ? 0 : -1;
			sh_ri = (i == xpts-1) ? 0 :  1;
			sh_do = (j == 0     ) ? 0 : -1;
			sh_up = (j == ypts-1) ? 0 :  1;

			p1 = p(i+sh_ri,j); p2 = p(i,j+sh_up); p3 = p(i+sh_le,j); p4 = p(i,j+sh_do);
			p5 = p(i+sh_ri,j+sh_up); p6 = p(i+sh_le,j+sh_up);
			p7 = p(i+sh_le,j+sh_do); p8 = p(i+sh_ri,j+sh_do);

			q1 = q(i+sh_ri,j); q2 = q(i,j+sh_up); q3 = q(i+sh_le,j); q4 = q(i,j+sh_do);
			q5 = q(i+sh_ri,j+sh_up); q6 = q(i+sh_le,j+sh_up);
			q7 = q(i+sh_le,j+sh_do); q8 = q(i+sh_ri,j+sh_do);
			
			JJ(i,j) = (p1-p3)*(q2-q4) - (p2-p4)*(q1-q3) \
					   + q2*(p5-p6) - q4*(p8-p7) - q1*(p5-p8) + q3*(p6-p7) \
					   + p1*(q5-q8) - p2*(q5-q6) - p3*(q6-q7) + p4*(q8-q7);

		});

	};



	int err;
	auto evolution = [&] (ModelData& k_vort, ModelData& evo_vort) {
		vort_tmp2.copyData(evo_vort);
		EnuData::each_index([&](xsize i, xsize j) {
			vort_tmp2(i,j) *= bc_filter(i,j);
		});

		// inverse vort to psi
		
		if((err = jacobi_relaxation(vort_tmp2, psi, tmp, dx, dy, 1, relative_err, max_iter)) != 0) {
			fprintf(stderr, "[Error code %d]: Relaxation cannot converge within [%d] iterations.\n", err, max_iter);
			exit(-1);
		}

		// calculate u,v
		// grady(psi, u, dy); u *= -1.0f;
		// gradx(psi, v, dx);

		// calculate grad vort
		// gradx(evo_vort, gradx_vort, dx);
		// grady(evo_vort, grady_vort, dy);

		// calculate lap
		//lap(evo_vort, lap_vort, dx, dy);

		ArakawaJacobian(evo_vort,psi,k_vort);
		k_vort /= (12*dx*dx);

//		EnuData::each_index([&](xsize i, xsize j){
//			k_vort(i,j) = -u(i,j) * gradx_vort(i,j) - v(i,j) * grady_vort(i,j) + nu_vort * lap_vort(i,j);	
//		});
	
	};

	auto calculator = [&] (ModelData& output, ModelData& input, ModelData& k, T dt) {
		EnuData:: each_index([&](xsize i, xsize j) {
			output(i,j) = input(i,j) + k(i,j) * dt;		
		});
	};


	stringstream vort_f;
	vort_f  << datafolder << "/vort_0.bin";
	output_data(vort , vort_f.str().c_str());

	time_t tstart, tend;
	tstart = time(0);
	for(xsize i = 1; i <= step; ++i) {
		cout << "Doing step " << i << ", time = " << dt * i << flush << ", ";
		evolution(k1_vort, vort);
			calculator(vort_tmp, vort, k1_vort, 0.5 * dt);
		evolution(k2_vort, vort_tmp);
			calculator(vort_tmp, vort, k2_vort, 0.5 * dt);
		evolution(k3_vort, vort_tmp);
			calculator(vort_tmp, vort, k3_vort,  dt);
		evolution(k4_vort, vort_tmp); 	 
		
		// calculate new vort
		EnuData::each_index([&](xsize i, xsize j){
			vort(i,j)  +=  dt / 6.0 * (k1_vort(i,j)  + 2*k2_vort(i,j)  + 2*k3_vort(i,j)  + k4_vort(i,j));
		});


		if(i % record_step == 0) {
			cout << "Output data timestep: " << i;
			vort_f.str("");
			vort_f  << datafolder <<"/vort_"  << i << ".bin";
			output_data(vort , vort_f.str().c_str());
		}

		cout << endl;
	}
	tend = time(0);
	cout << "It took " << difftime(tend, tstart) << " second(s)." << endl;
	return 0;
}


