#ifndef __AVERAGE_LIPHE__
#define __AVERAGE_LIPHE__

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
//	CostFunction *_out_cost;

public:
	AverageLiphe(int n) : _n(n), _avg(BinomialTournament<OutNumber>::add) {}
//	AverageLiphe(int n, CostFunction *c) : _n(n), _avg(BinomialTournament<OutNumber>::add), _cost(c) {}

	void add(const InNumber &x) {
		int a_i = random() % _n;
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

	void add_with_cost(const InNumber &x, const OutNumber &cost) {
		int a_i = random() % _n;
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

	OutNumber getAverage() const { return _avg.unite_all(); }
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
