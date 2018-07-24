#include "vector_operation.hpp"

#ifndef VECTOR_OPERATION_CPP
#define VECTOR_OPERATION_CPP
template<int PTS> void VectorOperation<PTS>::set(void *dat, void *val) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)dat)[i] = ((float *)val)[i];
	}
}

template<int PTS> void VectorOperation<PTS>::set(void *dat, float k) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)dat)[i] = k;
	}
}



template<int PTS> void VectorOperation<PTS>::add(void *out, void *in1, void *in2) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in1)[i] + ((float *)in2)[i];
	}
}

template<int PTS> void VectorOperation<PTS>::sub(void *out, void *in1, void *in2) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in1)[i] - ((float *)in2)[i];
	}
}

template<int PTS> void VectorOperation<PTS>::mul(void *out, void *in1, void *in2) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in1)[i] * ((float *)in2)[i];
	}
}

template<int PTS> void VectorOperation<PTS>::div(void *out, void *in1, void *in2) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in1)[i] / ((float *)in2)[i];
	}
}


// With scalar

template<int PTS> void VectorOperation<PTS>::add(void *out, void *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in)[i] + k;
	}
}
template<int PTS> void VectorOperation<PTS>::sub(void *out, void *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in)[i] - k;
	}
}
template<int PTS> void VectorOperation<PTS>::mul(void *out, void *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in)[i] * k;
	}
}
template<int PTS> void VectorOperation<PTS>::div(void *out, void *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		((float *)out)[i] = ((float *)in)[i] / k;
	}
}


#endif
