#ifndef ___ZP____
#define ___ZP____

#include <assert.h>
#include <algorithm>
#include <istream>
#include <ostream>
#include <vector>
#include <functional>

#include <stdio.h>

template<int SIMD_SIZE = 1>
class ZP {
private:
	static long long _prev_p;
	static long long _prev_r;
	long long _p;
	long long _r;
	long long _val[SIMD_SIZE];

	int _mul_depth;
	int _add_depth;

	long long find_inv(long long x, long long p, long long r) const;

	static long long power(long long base, long long exp) {
		long long ret = 1;
		for (int i = 0; i < exp; ++i)
			ret *= base;
		return ret;
	}

	void set_all(long long a) { for (int i = 0; i < SIMD_SIZE; ++i) _val[i] = mod(a); }
		
public:
	ZP() : _p(_prev_p), _r(_prev_r),  _mul_depth(0), _add_depth(0) {}
	ZP(long long v) : _p(_prev_p), _r(_prev_r),  _mul_depth(0), _add_depth(0) { set_all(v); }
	ZP(const std::vector<long int> &v) : _p(_prev_p), _r(_prev_r),  _mul_depth(0), _add_depth(0) { for (unsigned int i = 0; (i < v.size()) && (i < SIMD_SIZE); ++i) _val[i] = v[i]; }
	ZP(long long v, long long p) : _p(p), _r(_prev_r),  _mul_depth(0), _add_depth(0) { set_all(v); }
	ZP(long long v, long long p, long long r) : _p(p), _r(r),  _mul_depth(0), _add_depth(0) { set_all(v); }
	ZP(const ZP &zp) : _p(zp._p), _r(zp._r), _mul_depth(zp._mul_depth), _add_depth(zp._add_depth) {
		for (int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = zp._val[i];
	}

	int in_range(int a) const { while (a < 0) a += power(_p,_r); return a % power(_p,_r); }
	static int static_in_range(int a) { while (a < 0) a += power(_prev_p,_prev_r); return a % power(_prev_p,_prev_r); }
	int to_int() const { return _val[0] % power(_p,_r); }
	std::vector<long int> to_vector() const { return std::vector<long int>(_val, _val + SIMD_SIZE); }
	ZP clone() { return ZP(_val, _p, _r); }

	ZP from_int(int i) const { return ZP(i); }
	static ZP static_from_int(int i) { return ZP(i); }

	static unsigned int simd_factor() { return SIMD_SIZE; }
	static void set_global_p(long long p, long long r = 1) { _prev_p = p; _prev_r = r; }
	static int global_p() { return _prev_p; }
	static int get_global_ring_size() { return power(_prev_p, _prev_r); }
	void set_p(long long p, long long r = 1) { _p = _prev_p = p; _r = _prev_r = r; }
	long long p() const { return _p; }
	long long r() const { return _r; }
	long long get_ring_size() const { return power(p(), r()); }
	long long inv(long long x) const { assert(r() == 1); return (p() <= 3) ? x : power(x, p() - 2); }
	long long mod(long long x) const { return ((x % get_ring_size()) + get_ring_size()) % get_ring_size(); }

	int add_depth() const { return _add_depth; }
	int mul_depth() const { return _mul_depth; }

	ZP &operator=(const ZP &b) {
		_p = b._p;
		_r = b._r;
		_mul_depth = b._mul_depth;
		_add_depth = b._add_depth;
		for (int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = b._val[i];
		return *this;
	}

	void divide_by_p() {
		assert(_r > 1);
		for (int i = 0; i < SIMD_SIZE; ++i)
			_val[i] /= _p;
	}

	ZP rotate_left(int step) const {
		ZP ret(*this);

		for (int i = 0; i < SIMD_SIZE; ++i)
			ret._val[i] = _val[(i + step) % SIMD_SIZE];
		return ret;
	}

	ZP operator-() const {
		ZP zp(*this);
		for (int i = 0; i < SIMD_SIZE; ++i)
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
		for (int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = mod(_val[i] - z._val[i]);
		_mul_depth = std::max(_mul_depth, z._mul_depth);
		_add_depth = std::max(_add_depth, z._add_depth) + 1;
	}

	void operator+=(const ZP &z) {
		assert(_p == z._p);
		assert(_r == z._r);
		for (int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = mod(_val[i] + z._val[i]);
		_mul_depth = std::max(_mul_depth, z._mul_depth);
		_add_depth = std::max(_add_depth, z._add_depth) + 1;
	}

	void operator*=(const ZP &z) {
		assert(_p == z._p);
		assert(_r == z._r);
		for (int i = 0; i < SIMD_SIZE; ++i)
			_val[i] = mod(_val[i] * z._val[i]);
		_mul_depth = std::max(_mul_depth, z._mul_depth) + 1;
		_add_depth = std::max(_add_depth, z._add_depth);
	}

	template<class BITS>
	BITS to_digits() const {
		BITS ret;
		ret.set_bit_length(r());
		for (int i = 0; i < r(); ++i) {
			ZP<SIMD_SIZE> bit;
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

		for (int i = 0; i < SIMD_SIZE; ++i) {
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


	template<int A>
	friend std::ostream &operator<<(std::ostream &out, const ZP<A> &z);

	template<int A>
	friend std::istream &operator>>(std::istream &out, ZP<A> &z);
};

template<int SIMD_SIZE>
long long ZP<SIMD_SIZE>::_prev_p = 2;

template<int SIMD_SIZE>
long long ZP<SIMD_SIZE>::_prev_r = 1;


template<int SIMD_SIZE>
inline std::ostream &operator<<(std::ostream &out, const ZP<SIMD_SIZE> &z) {

	out
		<< z._p << " "
		<< z._r << " ";

	out << SIMD_SIZE << " ";
	for (int i = 0; i < SIMD_SIZE; ++i)
		out << z._val[i] << " ";

	out
		<< z._mul_depth << " "
		<< z._add_depth << " ";

	return out;
}

template<int SIMD_SIZE>
inline std::istream &operator>>(std::istream &in, ZP<SIMD_SIZE> &z) {
	in
		>> z._p
		>> z._r;

	int simd_size;
	in >> simd_size;
	assert(simd_size == SIMD_SIZE);
	for (int i = 0; i < SIMD_SIZE; ++i)
		in >> z._val[i];

	in
		>> z._mul_depth
		>> z._add_depth;

	return in;
}

#endif
