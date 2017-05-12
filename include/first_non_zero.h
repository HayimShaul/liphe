#ifndef ___FIRST_NON_ZERO__
#define ___FIRST_NON_ZERO__

#include <iostream>
#include <carray_iterator.h>
#include <binomial_tournament.h>

template<class Number, class IsZeroStrategy>
class BinarySearch {
public:
	// gets a binary array @array, 
//	static void searchFirst(const std::vector<Number> &array, IsZeroFunction *IsZero, std::vector<Number> &output) {
//		BinomialTournament>Number> bins;
//
//		for (int i = 0; i < output.size(); ++i)
//			output[i] = 0;
//
//		// No need to iterate over i == 0 because we are going to multiply everything by i==0
//		bins.add(array[0]);
//		for (int i = 1; i < size; ++i) {
//			prev = isZero(bin.unite_all(bin::mul));
//			int l = 0;
//			while ((1<<l) <= i) {
//				if (i & (1<<l))
//					output[l] += array[i] * prev * i;
//				++l;
//			}
//		}
//	}

	static void searchFirst(Number *array, int size, std::vector<Number> &output) {
		CArrayIterator<Number> begin(array, size, 0);
		CArrayIterator<Number> end(array, size, size);
		searchFirst(begin, end, output);
	}

	template<class Iterator>
	static void searchFirst(const Iterator &begin, const Iterator &end, std::vector<Number> &_output) {
		BinomialTournament<Number> bins(BinomialTournament<Number>::add);

		AddBinomialTournament<Number> *output = new AddBinomialTournament<Number>[_output.size()];

		// No need to iterate over i == 0 because we are going to multiply everything by i==0
		Iterator i = begin;
		bins.add_to_tournament(*i);
		++i;
		int count = 1;
		while (i != end) {
			BinomialTournament<Number> mul(BinomialTournament<Number>::mul);
			for (int level = 0; level < bins.max_level(); ++level) {
				if (!bins.is_slot_empty(level)) {
					// TODO replace with the new COMPARE standard in cmp.h
					mul.add_to_tournament(IsZeroStrategy::is(bins.number(level)));
				}
			}
			// TODO replace with the new COMPARE standard in cmp.h
			Number prev = IsZeroStrategy::is_not(mul.unite_all());

			for (int bit = 0; bit < _output.size(); ++bit) {
				if (count & (1 << bit)) {
					output[bit].add_to_tournament(prev * (*i));
				}
			}
			bins.add_to_tournament(*i);

			++count;
			++i;
		}
		for (int bit = 0; bit < _output.size(); ++bit) {
			// TODO replace with the new COMPARE standard in cmp.h
			_output[bit] = IsZeroStrategy::is_not(output[bit].unite_all());
		}
	}
};

#endif
