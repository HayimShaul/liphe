#ifndef __AVERAGE_LIPHE__
#define __AVERAGE_LIPHE__

#include <math.h>
#include <vector>
#include <random>

#include <print.h>
#include "binomial_tournament.h"



template<class OutNumber, class InNumber, class CompareIn, class Convert>
class AverageLiphe {
public:
//	typedef OutNumber OutCostFunction(const OutNumber &);
//	typedef InNumber OutCostFunction(const OutNumber &);

private:

	int _m;
	int _n;
	int _c;
	int _max_sample;

	std::function<int(int)> f_m;

	BinomialTournament<OutNumber> _avg;
	int _max_cost;
//	CostFunction *_out_cost;

	std::default_random_engine _generator;
	std::uniform_int_distribution<int> _distribution;

	int ceiling(float f) { return (int)(f+0.9999999); }

public:
//	AverageLiphe(int n) : _m(n), _n(n), _avg(BinomialTournament<OutNumber>::add), _c(1), _max_cost(1), _generator(clock()), _distribution(0, n) {}
	AverageLiphe(int n) : _m(n), _n(n), _c(1), _max_sample(0), _avg(BinomialTournament<OutNumber>::add), _max_cost(1), _generator(1), _distribution(0, n)  {}
//	AverageLiphe(int n, CostFunction *c) : _n(n), _avg(BinomialTournament<OutNumber>::add), _cost(c) {}

	void set_resample_constant(int r) { _c = r; }
	void set_m(int m) { _m = m; _distribution = std::uniform_int_distribution<int>(0, m); }
	void set_max_cost(int c) { _max_cost = c; }
	void set_max_sample(int c) { _max_sample = c; }
	void set_f_m(const std::function<int(int)> &f) { f_m = f; }

	void compute_resample_constant(float delta, float epsilon) {
		assert(_max_sample > 0);
		_c = 1;

		// make sure [cm]/sample < 1
		if ((bool)f_m) {
//			for (int m_i = 0; i <= _m; ++i)
//				_c = std::max(_c, ceiling( (float) max_sample / f_m(m_i) ) );
			// assume f_m is monotonously increasing
			_c = std::max(_c, ceiling( (float) _max_sample / f_m(_m) ) );
		} else {
			_c = std::max(_c, ceiling( (float) _max_sample / _m ) );
		}

//		// convert between Pr(...>delta) < eps    to   Pr(...>delta) < 2^-eps
//		epsilon = - log(epsilon) / log(2);
//		// make sure pr(..>delta) < eps
//		_c = std::max(_c, ceiling( (float) _max_cost * epsilon / (delta * delta * 2 * _n)));

		_distribution = std::uniform_int_distribution<int>(0, _c * _m); 

		std::cerr << "c = " << _c << std::endl;
	}

	void add(const InNumber &x) {
		for (int i = 0; i < _c; ++i) {
			int a_i = _distribution(_generator);
			if ((bool) f_m)
				a_i = f_m(a_i);

			// the case of a_i=0 is not interesting to compare to because when x_i = a_i = 0 it is not
			// going to contribute to the average anyway

			OutNumber addon = Convert::convert(CompareIn(x) > a_i);

			_avg.add_to_tournament( addon );
		}
	}

	void add_simd(const InNumber &x) {
		for (int i = 0; i < _c; ++i) {
			std::vector<int> a_i(x.simd_factor());

			for (unsigned int simd_i = 0; simd_i < a_i.size(); ++simd_i) {
				a_i[simd_i] = _distribution(_generator);
				if ((bool) f_m)
					a_i[simd_i] = f_m(a_i[simd_i]);
				if ((_max_sample > 0) && (a_i[simd_i] > _max_sample)) {
					a_i[simd_i] = _max_sample + 1;
				}
			}

			InNumber encAi(a_i);

			// the case of a_i=0 is not interesting to compare to because when x_i = a_i = 0 it is not
			// going to contribute to the average anyway

			OutNumber addon = Convert::convert(CompareIn(x) > encAi);

			_avg.add_to_tournament( addon );
		}
	}

	void add_with_cost(const InNumber &x, const OutNumber &cost) {
		for (int i = 0; i < _c; ++i) {
			int a_i = random() % (_m * _c);
			if ((bool) f_m)
				a_i = f_m(a_i);
			// the case of a_i=0 is not interesting to compare to because when x_i = a_i = 0 it is not
			// going to contribute to the average anyway

//std::cerr << "comparing   *i > a_1 = " << (*i) << " > " << a_i << std::endl;
//std::cerr << "result = " << (CompareIn(*i) > a_i) << std::endl;
			OutNumber addon = Convert::convert(CompareIn(x) > a_i);
//std::cerr << "addon = " << addon << std::endl;

//int a = addon.to_int();
//std::cerr << " adding  " << a << std::endl;

			_avg.add_to_tournament( addon * cost );
		}
	}

	OutNumber getAverage() const {
		if (_avg.is_empty())
			return OutNumber(0);
		return _avg.unite_all();
	}
};







// InNumber is usually the bit-wise reprsenation of xi
// OutNumber is a native number in Zp
// UnsignedWord is UnsignedWord<OutNumber>  normally that would be the same as InNuber
template<class UnsignedWord, class OutNumber, class InNumber, class CompareIn, class Convert>
class AverageLipheBits {
	int _n;
	int _m;
	int _max_cost;
	int _max_sample;
	float _delta;
	float _epsilon;
	std::vector< AverageLiphe<OutNumber, InNumber, CompareIn, Convert> * > _bits;

	void clean() {
		for (int i = 0; i != _bits.size(); ++i) {
			if (_bits[i] != NULL) {
				delete _bits[i];
				_bits[i] = NULL;
			}
		}
	}

	void apply_constants() {
		for (int i = 0; i < _bits.size(); ++i) {
			_bits[i]->set_max_cost(_max_cost);

//			_bits[i]->set_m(_m * (1 << i));
			_bits[i]->set_m(_m);

			if ((_delta > 0) && (_epsilon > 0) && (_max_sample > 0))
				_bits[i]->compute_resample_constant(_delta, _epsilon, _max_sample);
		}
	}

public:
	AverageLipheBits(int n) : _n(n), _m(n), _max_cost(1), _max_sample(0), _delta(0), _epsilon(0)  {}
	~AverageLipheBits() { clean(); }

	void set_m(int m) { _m = m; apply_constants(); }
	void set_max_cost(int m) { _max_cost = m; apply_constants(); }

	void compute_resample_constant(float delta, float epsilon, int max_sample) {
		_delta = delta;
		_epsilon = epsilon;
		_max_sample = max_sample;
		apply_constants();
	}

	
	void set_bit_number(int b) {
		clean();
		_bits.resize(b);
		for (int i = 0; i < b; ++i)
			_bits[i] = new AverageLiphe<OutNumber, InNumber, CompareIn, Convert>(_m);
		apply_constants();
	}

	void add(const InNumber &xi) {
		for (int i = 0; i < _bits.size(); ++i) {
			InNumber _xi = xi;
			_xi >>= i;
//std::cerr << "xi = ";
//std::cerr << _xi.to_int() << std::endl;
			_bits[i]->add(_xi);
		}
	}

	void add_with_cost(const InNumber &xi, const OutNumber &cost) {
//		for (int i = 0; i < _bits.size(); ++i)
//			_bits[i]->add_with_cost(xi, cost);
	}

	UnsignedWord get_bits() const {
		UnsignedWord ret;
		ret.set_bit_length(_bits.size());

		int i;
		for (i = 0; i < _bits.size(); ++i)
			ret[i] = _bits[i]->getAverage();

		return ret;
	}

	UnsignedWord getAverage() const {
		UnsignedWord ret;
		ret.set_bit_length(_bits.size());

		std::vector<OutNumber> bits(_bits.size());
		int i;
		for (i = 0; i < _bits.size(); ++i)
			bits[i] = _bits[i]->getAverage();

		i = _bits.size() - 1;
		ret[i] = bits[i];
		--i;
		while (i >= 0) {
			ret[i] = bits[i] - bits[i+1]*2;
			--i;
		}

		return ret;
	}
};








template<class OutNumber, class CompareIn, class Convert, class Iterator>
inline OutNumber average_by_probability(const Iterator &begin, const Iterator &end) {
	// assert all numbers have the same p
	// assert (n > p)
	
	BinomialTournament<OutNumber> average( BinomialTournament<OutNumber>::add );

	int n = end - begin;
	Iterator i = begin;

	int count = 0;
	while (i != end) {
		std::cerr << "Doing point " << count << "\r";
		int a_i = random() % n;
		// the case of a_i=0 is not interesting to compare to because when x_i = a_i = 0 it is not
		// going to contribute to the average anyway

//std::cerr << "comparing   *i > a_1 = " << (*i) << " > " << a_i << std::endl;
//std::cerr << "result = " << (CompareIn(*i) > a_i) << std::endl;
		OutNumber addon = Convert::convert(CompareIn(*i) > a_i);
//std::cerr << "addon = " << addon << std::endl;

//int a = addon.to_int();
//std::cerr << " adding  " << a << std::endl;

		average.add_to_tournament( addon );

		++i;
		++count;
	}
	std::cerr << std::endl;
	return average.unite_all();
}

//compute by \sum a_i   where a_i is 1 with probability x_i/n
template<class OutNumber, class InNumber, class CompareIn>
inline OutNumber average_by_probability(InNumber *array, int size) {
	CArrayIterator<InNumber> begin(array, size, 0);
	CArrayIterator<InNumber> end(array, size, size);

	return average_by_probability<OutNumber, CompareIn>(begin, end);
}

#endif
