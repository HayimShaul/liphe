#ifndef __AVERAGE_LIPHE__
#define __AVERAGE_LIPHE__

#include <math.h>
#include <vector>

#include <print.h>
#include "binomial_tournament.h"



template<class OutNumber, class InNumber, class CompareIn, class Convert>
class AverageLiphe {
public:
//	typedef OutNumber OutCostFunction(const OutNumber &);
//	typedef InNumber OutCostFunction(const OutNumber &);

private:
	int _n;
	int _count;
	BinomialTournament<OutNumber> _avg;
	int _resample;
	int _n_factor;
	int _max_cost;
//	CostFunction *_out_cost;

	int ceiling(float f) { return (int)(f+0.9999999); }

public:
	AverageLiphe(int n) : _n(n), _avg(BinomialTournament<OutNumber>::add), _resample(1), _n_factor(1), _max_cost(1) {}
//	AverageLiphe(int n, CostFunction *c) : _n(n), _avg(BinomialTournament<OutNumber>::add), _cost(c) {}

	void set_resample_constant(int r) { _resample = r; }
	void set_n_factor(int r) { _n_factor = r; }
	void set_max_cost(int c) { _max_cost = c; }

	void compute_resample_constant(float delta, float epsilon, int max_sample) {
		_resample = 1;

		int c;
		// make sure [cm]/sample < 1
		_resample = std::max(_resample, ceiling( (float) max_sample / (_n * _n_factor) ) );

		// convert between Pr(...>delta) < eps    to   Pr(...>delta) < 2^-eps
		epsilon = - log(epsilon) / log(2);
		// make sure pr(..>delta) < eps
		_resample = std::max(_resample, ceiling( (float) _max_cost * epsilon / (2 * _n * delta * delta)));

		std::cerr << "resample = " << _resample << std::endl;
	}

	void add(const InNumber &x) {
		for (int i = 0; i < _resample; ++i) {
			int a_i = random() % (_n * _n_factor * _resample);
			// the case of a_i=0 is not interesting to compare to because when x_i = a_i = 0 it is not
			// going to contribute to the average anyway

//std::cerr << "comparing   *i > a_1 = " << (*i) << " > " << a_i << std::endl;
//std::cerr << "result = " << (CompareIn(*i) > a_i) << std::endl;
			OutNumber addon = Convert::convert(CompareIn(x) > a_i);
//std::cerr << "addon = " << addon << std::endl;

//int a = addon.to_int();
//std::cerr << " adding  " << a << std::endl;

//		assert(_cost == NULL);
			_avg.add_to_tournament( addon );
		}
	}

	void add_with_cost(const InNumber &x, const OutNumber &cost) {
		for (int i = 0; i < _resample; ++i) {
			int a_i = random() % (_n * _resample);
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

	OutNumber getAverage() const { return _avg.unite_all(); }
};







// InNumber is usually the bit-wise reprsenation of xi
// OutNumber is a native number in Zp
// UnsignedWord is UnsignedWord<OutNumber>  normally that would be the same as InNuber
template<class UnsignedWord, class OutNumber, class InNumber, class CompareIn, class Convert>
class AverageLipheBits {
	int _n;
	int _max_cost;
	int _max_sample;
	int _n_factor;
	float _delta;
	float _epsilon;
	std::vector< AverageLiphe<OutNumber, InNumber, CompareIn, Convert> * > _bits;
	std::vector<int> _resample;

	void clean() {
		for (int i = 0; i != _bits.size(); ++i) {
			if (_bits[i] != NULL) {
				delete _bits[i];
				_bits[i] = NULL;
			}
		}
	}

	void apply_constants() {
		assert((_resample.size() == _bits.size()) || (_resample.size() == 0));
		for (int i = 0; i < _bits.size(); ++i) {
			if (_resample.size() == 0)
				_bits[i]->set_resample_constant(1);
			else
				_bits[i]->set_resample_constant(_resample[i]);

			_bits[i]->set_max_cost(_max_cost);

			_bits[i]->set_n_factor(_n_factor * (1 << i));

			if ((_delta > 0) && (_epsilon > 0) && (_max_sample > 0))
				_bits[i]->compute_resample_constant(_delta, _epsilon, _max_sample);
		}
	}

public:
	AverageLipheBits(int n) : _n(n), _max_cost(1), _delta(0), _epsilon(0), _n_factor(1)  {}
	~AverageLipheBits() { clean(); }

	void set_n_factor(int m) { _n_factor = m; apply_constants(); }
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
		for (int i = 0; i < b; ++i) {
			_bits[i] = new AverageLiphe<OutNumber, InNumber, CompareIn, Convert>(_n);
			_bits[i]->set_n_factor(1 << i);
		}
		apply_constants();
	}

	void set_resample_constant(std::vector<int> c) {
		_resample = c;
		apply_constants();
	}

	void add(const InNumber &xi) {
		for (int i = 0; i < _bits.size(); ++i)
			_bits[i]->add(xi);
	}

	void add_with_cost(const InNumber &xi, const OutNumber &cost) {
		for (int i = 0; i < _bits.size(); ++i)
			_bits[i]->add_with_cost(xi, cost);
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
