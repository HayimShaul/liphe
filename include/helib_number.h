#ifndef ___HELIB_NUMBER___
#define ___HELIB_NUMBER___

#include <assert.h>
#include <helib_keys.h>

class HelibNumber {
private:
	static HelibKeys *_prev_keys;
	HelibKeys *_keys;
	Ctxt _val;

	int _mul_depth;
	int _add_depth;

	long long power(long long base, long long exp) const {
		long long ret = 1;
		for (int i = 0; i < exp; ++i)
			ret *= base;
		return ret;
	}

	HelibNumber mult_by_recursive_adding(const HelibNumber &x, int e) {
		assert(e >= 0);
		if (e == 0)
			return HelibNumber(0);

		if (e == 2) {
			HelibNumber y = x + x;
			return y;
		}
		if (e == 1)
			return x;

		if (e & 1) {
			return x + mult_by_recursive_adding(x, e-1);
		} else {
			return mult_by_recursive_adding(x + x, e/2);
		}
	}



public:
	HelibNumber() : _keys(_prev_keys), _val(_prev_keys->publicKey()), _mul_depth(0), _add_depth(0) {}
	HelibNumber(long long v) : _keys(_prev_keys), _val(_keys->publicKey()), _mul_depth(0), _add_depth(0) { _keys->encrypt(_val, v); }
	HelibNumber(const HelibNumber &n) : _keys(n._keys), _val(n._val), _mul_depth(n._mul_depth), _add_depth(n._add_depth) {}

	static void set_global_keys(HelibKeys *k) { _prev_keys = k; }
	int in_range(int a) const { while (a < 0) a += _keys->p(); return a % _keys->p(); }
	static int static_in_range(int a) { while (a < 0) a += _prev_keys->p(); return a % _prev_keys->p(); }
	int to_int() const { return _keys->decrypt(_val); }
	HelibNumber from_int(int i) const {
			HelibNumber ret;
			_keys->encrypt(ret._val, i);
			return ret;
		}

	static HelibNumber static_from_int(int i) {
			HelibNumber ret;
			_prev_keys->encrypt(ret._val, i);
			return ret;
		}

	int get_ring_size() const { return power(_keys->p(), _keys->r()); }
	int p() const { return _keys->p(); }
	int r() const { return _keys->r(); }
	void assert_co_prime(int) const {}
//	HelibNumber clone() { return HelibNumber(*this); }

	HelibNumber &operator=(const HelibNumber &b) {
		_keys = b._keys;
		_val = b._val;
		_mul_depth = b._mul_depth;
		_add_depth = b._add_depth;
		return *this;
	}

	const Ctxt &v() const { return _val; }

	int add_depth() const { return _add_depth; }
	int mul_depth() const { return _mul_depth; }

	void shift_right() { _val.divideByP(); }

	void negate() { _val.negate(); }




	HelibNumber operator-() const { HelibNumber zp(*this); zp.negate(); return zp; }
	HelibNumber operator-(const HelibNumber &z) const { HelibNumber zp(*this); zp -= z; return zp; }
	HelibNumber operator+(const HelibNumber &z) const { HelibNumber zp(*this); zp += z; return zp; }
	HelibNumber operator*(const HelibNumber &z) const { HelibNumber zp(*this); zp *= z; return zp; }

	void operator-=(const HelibNumber &z) {
		assert(_keys == z._keys);
//std::cerr << "before -=\n";
//std::cerr << "  level=" << _val.findBaseLevel() << ", log(noise/modulus)~" << _val.log_of_ratio() << endl;
		_val -= z._val;
//std::cerr << "after -=\n";
//std::cerr << "  level=" << _val.findBaseLevel() << ", log(noise/modulus)~" << _val.log_of_ratio() << endl;
		_mul_depth = max(_mul_depth, z._mul_depth);
		_add_depth = max(_add_depth, z._add_depth) + 1;
	}

	void operator+=(const HelibNumber &z) {
		assert(_keys == z._keys);
		_val += z._val;
		_mul_depth = max(_mul_depth, z._mul_depth);
		_add_depth = max(_add_depth, z._add_depth) + 1;
	}

	void operator*=(const HelibNumber &z) {
		assert(_keys == z._keys);
		_val.multiplyBy(z._val);
		_mul_depth = max(_mul_depth, z._mul_depth) + 1;
		_add_depth = max(_add_depth, z._add_depth);
	}






	HelibNumber operator-(int z) const { HelibNumber zp(*this); zp -= z; return zp; }
	HelibNumber operator+(int z) const { HelibNumber zp(*this); zp += z; return zp; }
	HelibNumber operator*(int z) const { HelibNumber zp(*this); zp *= z; return zp; }

	void operator-=(int z) { operator-=(HelibNumber(z)); }
	void operator+=(int z) { operator+=(HelibNumber(z)); }
	void operator*=(int z) { *this = mult_by_recursive_adding(*this, z); }


	void reduceNoiseLevel() {
		_val.modDownToLevel(_val.findBaseLevel());
	}
};

#endif
