#ifndef ___ZP____
#define ___ZP____

#include <assert.h>
#include <algorithm>
#include <istream>
#include <ostream>
#include <vector>
#include <functional>

#include <stdio.h>

class ZP {
private:
	static unsigned int SIMD_SIZE;

	static long _prev_p;
	static long _prev_r;

	static long getPrevP() { return _prev_p; }
	static long getPrevR() { return _prev_r; }

	static std::function<long (void)> _getP;
	static std::function<long (void)> _getR;

	long _p;
	long _r;
	long *_val;

	int _mul_depth;
	int _add_depth;

	long find_inv(long x, long p, long r) const;

	static long power(long base, long exp) {
		long ret = 1;
		for (int i = 0; i < exp; ++i)
			ret *= base;
		return ret;
	}

	void set_all(long a) { for (unsigned int i = 0; i < SIMD_SIZE; ++i) _val[i] = in_range(a); }
	void set_all(long *a) { for (unsigned int i = 0; i < SIMD_SIZE; ++i) _val[i] = in_range(a[i]); }
		
public:
	ZP() : _p(_getP()), _r(_getR()),  _mul_depth(0), _add_depth(0) { _val = new long[SIMD_SIZE]; }
	ZP(long v) : _p(_getP()), _r(_getR()), _mul_depth(0), _add_depth(0) { _val = new long[SIMD_SIZE]; set_all(v); }
	ZP(const std::vector<long> &v) : _p(_getP()), _r(_getR()),  _mul_depth(0), _add_depth(0) { _val = new long[SIMD_SIZE]; for (unsigned int i = 0; (i < v.size()) && (i < SIMD_SIZE); ++i) _val[i] = in_range(v[i]); }
	ZP(long v, long p) : _p(p), _r(_getR()),  _mul_depth(0), _add_depth(0) { _val = new long[SIMD_SIZE]; set_all(v); }
	ZP(long v, long p, long r) : _p(p), _r(r),  _mul_depth(0), _add_depth(0) { _val = new long[SIMD_SIZE]; set_all(v); }
	ZP(long *v, long p) : _p(p), _r(_getR()),  _mul_depth(0), _add_depth(0) { _val = new long[SIMD_SIZE]; set_all(v); }
	ZP(long *v, long p, long r) : _p(p), _r(r),  _mul_depth(0), _add_depth(0) { _val = new long[SIMD_SIZE]; set_all(v); }
	ZP(const ZP &zp) : _p(zp._p), _r(zp._r), _mul_depth(zp._mul_depth), _add_depth(zp._add_depth) {
		_val = new long[SIMD_SIZE];
 		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = zp._val[i];
	}

	~ZP() { delete[] _val; }

	int in_range(int a) const { while (a < 0) a += power(_p,_r); return a % power(_p,_r); }
	static int static_in_range(int a) { while (a < 0) a += power(_getP(),_getR()); return a % power(_getP(),_getR()); }
	int to_int() const { return _val[0] % power(_p,_r); }
	std::vector<long> to_vector() const { return std::vector<long>(_val, _val + SIMD_SIZE); }
	ZP clone() { return ZP(_val, _p, _r); }

	void from_int(long i) { _val[0] = in_range(i); }
	static ZP static_from_int(long i) { return ZP(i); }
	void from_vector(const std::vector<long> &in) {
		for (unsigned int i = 0; i < std::min((unsigned long)SIMD_SIZE, in.size()); ++i)
			_val[i] = in_range(in[i]);
	}

	static unsigned int simd_factor() { return SIMD_SIZE; }
	static void set_global_simd_factor(int a) { SIMD_SIZE = a; }
	static void set_global_p(long p, long r = 1) { _prev_p = p; _prev_r = r; }
	static int global_p() { return _getP(); }
	static int get_global_ring_size() { return power(_getP(), _getR()); }
	void set_p(long p, long r = 1) { _p = _prev_p = p; _r = _prev_r = r; }
	long p() const { return _p; }
	long r() const { return _r; }
	long get_ring_size() const { return power(p(), r()); }
	long inv(long x) const { assert(r() == 1); return (p() <= 3) ? x : power(x, p() - 2); }
	long mod(long x) const { return ((x % get_ring_size()) + get_ring_size()) % get_ring_size(); }

	int add_depth() const { return _add_depth; }
	int mul_depth() const { return _mul_depth; }

	ZP &operator=(const ZP &b) {
		_p = b._p;
		_r = b._r;
		_mul_depth = b._mul_depth;
		_add_depth = b._add_depth;
		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = b._val[i];
		return *this;
	}

	void divide_by_p() {
		assert(_r > 1);
		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			_val[i] /= _p;
	}

	ZP shift_left(int step) const {
		ZP ret(*this);

		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			if (i + step < SIMD_SIZE)
				ret._val[i] = _val[i + step];
		return ret;
	}

	ZP shift_right(int step) const {
		ZP ret(*this);

		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			if (i - step >= 0)
				ret._val[i] = _val[i - step];
		return ret;
	}

	ZP rotate_left(int step) const {
		ZP ret(*this);

		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			ret._val[i] = _val[(i + step) % SIMD_SIZE];
		return ret;
	}

	ZP operator-() const {
		ZP zp(*this);
		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			zp._val[i] = -zp._val[i];
		return zp;
	}

	ZP operator-(const ZP &z) const { ZP zp(*this); zp -= z; return zp; }
	ZP operator+(const ZP &z) const { ZP zp(*this); zp += z; return zp; }
	ZP operator*(const ZP &z) const { ZP zp(*this); zp *= z; return zp; }
	ZP operator!() const { ZP zp(*this); zp = ZP(1) - zp; return zp; }

	void operator-=(const ZP &z) {
		assert(_p == z._p);
		assert(_r == z._r);
		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = mod(_val[i] - z._val[i]);
		_mul_depth = std::max(_mul_depth, z._mul_depth);
		_add_depth = std::max(_add_depth, z._add_depth) + 1;
	}

	void operator+=(const ZP &z) {
		assert(_p == z._p);
		assert(_r == z._r);
		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = mod(_val[i] + z._val[i]);
		_mul_depth = std::max(_mul_depth, z._mul_depth);
		_add_depth = std::max(_add_depth, z._add_depth) + 1;
	}

	void operator*=(const ZP &z) {
		assert(_p == z._p);
		assert(_r == z._r);
		for (unsigned int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = mod(_val[i] * z._val[i]);
		_mul_depth = std::max(_mul_depth, z._mul_depth) + 1;
		_add_depth = std::max(_add_depth, z._add_depth);
	}

	template<class BITS>
	BITS to_digits() const {
		BITS ret;
		ret.set_bit_length(r());
		for (int i = 0; i < r(); ++i) {
			ZP bit;
			for (int s = 0; s < SIMD_SIZE; ++s) {
				bit._val[s] = (_val[s] / power(p(), i)) % p();
			}
			ret.set_bit(i, bit);
		}
		return ret;
	}

	void assert_co_prime(int a) const {
		assert(a != 0);
		assert(a != 1);

		for (unsigned int i = 0; i < SIMD_SIZE; ++i) {
			int b = _val[i];

			if (b > 1) {

				while (((b % a) != 0) && ((b % a) != 0)) {
					if (a > b)
						a -= b;
					else
						b -= a;
				}
				int gcd = (a < b) ? a : b;
				assert (gcd == 1);
			}
		}
	}

//	bool operator>(const ZP &z) const { assert(_p == z._p); return _val[0] > z._val[0]; }
//	bool operator<(const ZP &z) const { assert(_p == z._p); return _val[0] < z._val[0]; }
//	bool operator!=(const ZP &z) const { assert(_p == z._p); return _val[0] != z._val[0]; }

	void reduceNoiseLevel() {}

	friend std::ostream &operator<<(std::ostream &out, const ZP &z);
	friend std::istream &operator>>(std::istream &out, ZP &z);
};


inline std::ostream &operator<<(std::ostream &out, const ZP &z) {

	out
		<< z._p << " "
		<< z._r << " ";

	out << z.simd_factor() << " ";
	for (unsigned int i = 0; i < z.simd_factor(); ++i)
		out << z._val[i] << " ";

	out
		<< z._mul_depth << " "
		<< z._add_depth << " ";

	return out;
}

inline std::istream &operator>>(std::istream &in, ZP &z) {
	in
		>> z._p
		>> z._r;

	unsigned int simd_size;
	in >> simd_size;
	assert(simd_size == z.simd_factor());
	for (unsigned int i = 0; i < z.simd_factor(); ++i)
		in >> z._val[i];

	in
		>> z._mul_depth
		>> z._add_depth;

	return in;
}

#endif
