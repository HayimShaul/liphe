#ifndef __TEST_FRAMEWORK__
#define __TEST_FRAMEWORK__

#include <iostream>

#include <stdlib.h>
#include <time.h>

#define skipDoTest(name, func, Number, data, iter, seed)                                                  \
	std::cout << "Skipping test " << name << std::endl;


#define doTest(name, func, Number, data, iter, seed)                                                      \
	do {                                                                                                 \
		srand(seed);                                                                                     \
		int i = 0;                                                                                       \
		bool ok = true;                                                                                  \
		clock_t start = clock();																		 \
		while ((i < iter) && (ok)) {                                                                     \
			ok = func<Number>(data);                                                                      \
			++i;                                                                                         \
		}                                                                                                \
		std::cout << "Test " << name << "     ";                                                         \
		if (ok) {                                                                                        \
			std::cout << "OK (" << (clock() - start) << " clocks)" << std::endl;                         \
		} else {                                                                                         \
			std::cout << "FAILED    (seed=" << seed << "  iter=" << i << ")" << std::endl;               \
			exit(0);                                                                                     \
		}                                                                                                \
	} while(false)


// testing arithmetic functionality

template<class Number>
inline bool test_add(void *data) {
	int a = Number::static_in_range(rand());
	int b = Number::static_in_range(rand());

	Number A(a);
	Number B(b);

	Number C = A + B;

	int c = C.to_int();

	return (c == Number::static_in_range(a+b));
}

template<class Number>
inline bool test_sub(void *data) {
	int a = Number::static_in_range(rand());
	int b = Number::static_in_range(rand());

	Number A(a);
	Number B(b);

	Number C = A - B;

	int c = C.to_int();

	return (c == Number::static_in_range(a-b));
}

template<class Number>
inline bool test_mul(void *data) {
	int a = Number::static_in_range(rand());
	int b = Number::static_in_range(rand());

	Number A(a);
	Number B(b);

	Number C = A * B;

	int c = C.to_int();

	return (c == Number::static_in_range(a*b));
}

template<class Number>
inline bool test_bitwise_less_than(void *data) {
	int a = Number::static_in_range(rand());
	int b = Number::static_in_range(rand());

	Number A(a);

	Number C = A < b;

	int c = C.to_int();

	return (c == (a < b));
}

template<class Number>
inline bool test_bitwise_more_than(void *data) {
	int a = Number::static_in_range(rand());
	int b = Number::static_in_range(rand());

	Number A(a);

	Number C = A > b;

	int c = C.to_int();

	return (c == (a > b));
}

template<class Number>
inline bool test_euler_eq(void *data) {
	int a = Number::static_in_range(rand());
	int b = Number::static_in_range(rand());

	Number A(a);
	Number B(b);

	Number C = eq_euler(A, B);

	int c = C.to_int();

	return (c == ((a == b) ? 1 : 0));
}

template<class Number>
inline bool test_binomial_tournament(void *data) {
	int size = rand() % 1000;

	int a;
	Number A;

	int sum = 0;
	BinomialTournament<Number> binTour(BinomialTournament<Number>::add);

	for (int i = 0; i < size; ++i) {
		a = Number::static_in_range(rand());
		A = Number(a);

		sum += a;
		binTour.add_to_tournament(A);
	}

	Number C = binTour.unite_all();
	int c = C.to_int();

	return (c == Number::static_in_range(sum));
}

template<class Number>
inline bool test_find_first_in_array(void *data) {
	int log_size = 3;
	int size = rand() % (1 << log_size);
	size = 1 << log_size;

	Number *array = new Number[size];
	int *int_array = new int[size];
	std::vector<Number> output;
	output.resize(log_size);

int_array[0] = 0;
int_array[1] = 0;
int_array[2] = 0;
int_array[3] = 1;
int_array[4] = 0;
int_array[5] = 1;
int_array[6] = 0;
int_array[7] = 1;

	for (int i = 0; i < size; ++i) {
//		int_array[i] = rand() & 1;
		array[i] = Number(int_array[i]);
	}

	BinarySearch<Number, IsZeroEulerStrategy<Number> >::searchFirst(array, 5, output);

	int class_place = 0;
	for (int i = 0; i < log_size; ++i) {
		class_place += (1 << i) * output[i].to_int();
	}

	int real_place = 0;
	while ((real_place < size) && (int_array[real_place] == 0))
		++real_place;

	delete[] array;
	delete[] int_array;

	return (real_place == class_place);
}

//inline bool test_find_first_in_array(void *data) {
//	int n = 10;
//
//	HelibNumber *numbers = new HelibNumber[n];
//	int *_numbers = new int[n];
//	for (int i = 0; i < n; ++i) {
//		numbers[i] _numbers[i] = i;
//	}
//
//	HelibNumber res;
//	int res_int;
//
//	res = average_by_probability(numbers);
//	res_int = res.to_int();
//	std::cout << "average by probability:" << res_int << std::endl;
//	
//	delete[] numbers;
//	delete[] _nmbers
//}


#endif
