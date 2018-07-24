#ifndef VECTOR_OPERATION_HPP
#define VECTOR_OPERATION_HPP

template<int PTS> class VectorOperation {
private:
public:
	//VectorOperation();
	//~VectorOperation();
	void set(void *dat, void *val);
	void set(void *dat, float k);
	
	void add(void *out, void *in1, void *in2);
	void sub(void *out, void *in1, void *in2);
	void mul(void *out, void *in1, void *in2);
	void div(void *out, void *in1, void *in2);
	
	void add(void *out, void *in, float k);
	void sub(void *out, void *in, float k);
	void mul(void *out, void *in, float k);
	void div(void *out, void *in, float k);

	inline void iadd(void *out, void *in) { this->add(out, out, in); }
	inline void isub(void *out, void *in) { this->sub(out, out, in); }
	inline void imul(void *out, void *in) { this->mul(out, out, in); }
	inline void idiv(void *out, void *in) { this->div(out, out, in); }
	
	inline void iadd(void *out, float k) { this->add(out, out, k); }
	inline void isub(void *out, float k) { this->sub(out, out, k); }
	inline void imul(void *out, float k) { this->mul(out, out, k); }
	inline void idiv(void *out, float k) { this->div(out, out, k); }
	
};


#endif
