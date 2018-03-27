#ifndef ___TIME_MEASUREMENTS___
#define ___TIME_MEASUREMENTS___

#include <sys/time.h>

class TakeTimes {
	time_t _start;
	time_t _total;

	time_t _real_start;
	time_t _real_total;

public:
	TakeTimes() : _start(0), _total(0), _real_start(0), _real_total(0) {}

	void start() {
		assert(_start == 0);
		_start = clock();

		struct timeval tv;
		gettimeofday(&tv, NULL);

		assert(_real_start == 0);
		_real_start = tv.tv_sec;
	}

	std::string stats(const std::string &label) {
		std::stringstream str;
		str << label << " took " << _total << " cpu and " << _real_total << " time";
		if (_real_total > 0)
			str << ", which is a parallelization factor of " << (0.01*((int) (_total / 10000) / _real_total));
		str << std::endl;
		return str.str();
	}

	std::string end(const std::string &label) {
		time_t add = clock() - _start;

		_total += add;
		_start = 0;

		struct timeval tv;
		gettimeofday(&tv, NULL);

		time_t real_add = tv.tv_sec - _real_start;

		_real_total += real_add;
		_real_start = 0;

		std::stringstream str;
		str << label << " took " << (add / 1000000) << " cpu and " << real_add << " time";
		if (real_add > 0)
			str << ", which is a parallelization factor of " << (0.01*((int) (add / 10000) / real_add));
		str << std::endl;
		return str.str();
	}
};

class AutoTakeTimes {
private:
	std::string _s;
	TakeTimes _t;
public:
	AutoTakeTimes(const std::string &s) : _s(s) { _t.start(); }
	~AutoTakeTimes() { std::cout << _t.end(_s); }
};

#endif
