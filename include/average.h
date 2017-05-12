#ifndef __AVERAGE_LIPHE__
#define __AVERAGE_LIPHE__

#include "binomial_tournament.h"

template<class CLASS, class COMPARE, class Iterator>
inline CLASS average_by_probability(const Iterator &begin, const Iterator &end) {
	// assert all numbers have the same p
	// assert (n > p)
	
	BinomialTournament<CLASS> average( BinomialTournament<CLASS>::add );

	int n = end - begin;
	Iterator i = begin;

	while (i != end) {
		int a_i = random() % n;
		// the case of a_i=0 is not interesting to compare to because when x_i = a_i = 0 it is not
		// going to contribute to the average anyway

		CLASS addon = COMPARE(*i) > a_i;
		average.add_to_tournament( addon );

		++i;
	}
	return average.unite_all();
}

//compute by \sum a_i   where a_i is 1 with probability x_i/n
template<class CLASS, class COMPARE>
inline CLASS average_by_probability(CLASS *array, int size) {
	CArrayIterator<CLASS> begin(array, size, 0);
	CArrayIterator<CLASS> end(array, size, size);

	return average_by_probability<CLASS, COMPARE>(begin, end);
}


#endif
