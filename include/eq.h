#ifndef __EQ__
#define __EQ__

#include <iostream>


template<class CLASS>
CLASS power(const CLASS &x, int e) {
	assert(e > 0);
	if (e == 2) {
		CLASS y = x * x;
		return y;
	}
	if (e == 1)
		return x;

	if (e & 1) {
		return x*power(x, e-1);
	} else {
		return power(x*x, e/2);
	}
}

class Factorization {
	enum { MAX_PRIMES = 1024 };

	int primes[MAX_PRIMES];
	int multiplicity[MAX_PRIMES];
	int prime_no;

	int find_prime(int p) {
		int i;
		for (i = 0; i < prime_no; ++i)
			if (primes[i] == p)
				return i;
		
		++prime_no;
		assert(prime_no < MAX_PRIMES);
		multiplicity[i] = 0;
		primes[i] = p;
		return i;
	}

	void add_factor(int p) { ++multiplicity[find_prime(p)]; }
public:
	Factorization() : prime_no(0) {}
	Factorization(int p) : prime_no(0) { factor(p); }

	void factor(int x) {
		int primes[] = {
				2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
				131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
				269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
				421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
				587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
				743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
				919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997,
				1000003
   				,0 };

		int p_i = 0;
		while (primes[p_i] != 0) {
			while ((x % primes[p_i]) == 0) {
				add_factor(primes[p_i]);
				x /= primes[p_i];
			}
			++p_i;
		}
		assert(x == 1);
	}

	int factors() const { return prime_no; }
	int get_prime(int i) const { return primes[i]; }
	int get_multiplicity(int i) const { return multiplicity[i]; }
};

inline int phi(int p) {
	int phi = 1;

	Factorization f(p);
	for (int i = 0; i < f.factors(); ++i) {
		if (f.get_multiplicity(i) == 1) {
			 phi *= f.get_prime(i) - 1;
		} else {
			 phi *= power(f.get_prime(i), f.get_multiplicity(i)) - power(f.get_prime(i), f.get_multiplicity(i) - 1);
		}
	}

	return phi;
}

template<class CLASS>
CLASS isNonZero_Euler(const CLASS &r) {
	int p = r.get_ring_size();
	r.assert_co_prime(p);
	return power(r, phi(p));
}

template<class CLASS>
class IsNonZeroEulerStrategy {
public:
	static CLASS is(const CLASS &r) { return isNonZero_Euler(r); }
	static CLASS is_not(const CLASS &r) { return isZero_Euler(r); }
};


template<class CLASS>
CLASS isZero_Euler(const CLASS &r) { return r.from_int(1) - isNonZero_Euler(r); }

template<class CLASS>
class IsZeroEulerStrategy {
public:
	static CLASS is(const CLASS &r) { return isZero_Euler(r); }
	static CLASS is_not(const CLASS &r) { return isNonZero_Euler(r); }
};






template<class CLASS>
CLASS neq_euler(const CLASS &a, const CLASS &b) {
	assert(a.get_ring_size() == b.get_ring_size());
	CLASS r = a - b;

	return isNonZero_Euler(r);
}

template<class CLASS>
CLASS eq_euler(const CLASS &a, const CLASS &b) { return a.from_int(1) - neq_euler(a, b); }

#endif

