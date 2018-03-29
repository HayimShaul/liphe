#include <zp.h>
#include <unsigned_word.h>
#include <eq.h>
#include <binomial_tournament.h>
#include <first_non_zero.h>

#include "test_framework.h"

typedef ZP<127> MyZP;

typedef UnsignedWord<5, ZP<1> > MyUnsignedWord;

int main(int, char**) {
	MyZP::set_global_p(101);

	skipDoTest("operator + <ZP<127> >", test_add, MyZP, NULL, 1, -1);
	skipDoTest("operator - <ZP<127> >", test_sub, MyZP, NULL, 1, -1);
	skipDoTest("operator * <ZP<127> >", test_mul, MyZP, NULL, 1, -1);
	skipDoTest("euler_eq <ZP<127> >", test_euler_eq, MyZP, NULL, 1, -1);
	skipDoTest("BinomialTournament <ZP<127> >", test_binomial_tournament, MyZP, NULL, 1, -1);
	skipDoTest("FindFirstNonZero <ZP<127> >", test_find_first_in_array, MyZP, NULL, 1, -1);

	skipDoTest("operator < <UnsignedWord < ZP<127> >", test_bitwise_less_than, MyUnsignedWord, NULL, 100000, -1);
	doTest("operator > <UnsignedWord < ZP<127> >", test_bitwise_more_than, MyUnsignedWord, NULL, 100000, -1);

	MyZP::set_global_p(2);
	skipDoTest("operator + <UnsignedWord < ZP<127> >", test_add, MyUnsignedWord, NULL, 1, -1);
	skipDoTest("operator - <UnsignedWord < ZP<127> >", test_sub, MyUnsignedWord, NULL, 1, -1);
	skipDoTest("operator * <UnsignedWord < ZP<127> >", test_mul, MyUnsignedWord, NULL, 1, -1);

	MyZP::set_global_p(101);
	doTest("Polynomial < ZP<127> >", test_polynomial, MyZP, NULL, 1, -1);
}

