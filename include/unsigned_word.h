#ifndef ___UNSIGNED_WORD___
#define ___UNSIGNED_WORD___

#include <iostream>
#include <assert.h>
#include <bitset>
#include <vector>
#include <functional>

#include "binomial_tournament.h"

// Author: Hayim Shaul 2016.

// UnsignedWord.
// implements a positive number by its bit representation. Bit can be either ZP (with p=2) or HELibNumber (p = 2)
// UnsignedWords are implemented to have a variable number of bits, with a maximum defined to be MAX_BIT_NUM
//
// Note: UnsignedWord can be only be positive, negative numbers can be implemented with 2-complement but require
// a fixed number of bit representation, or some other TIHKUMATION like keeping the sign in a separate bit
// and for every + operation compute both +/- and choose by multiplying with sign bit

template<int MAX_BIT_NUM, class Bit>
class UnsignedWord {
private:
	std::vector<Bit *> _bits;

	Bit biggerBiggerEqual(const UnsignedWord<MAX_BIT_NUM, Bit> &w, bool isBiggerEqual) const;
public:
	static unsigned int simd_factor() { return Bit::simd_factor(); }

	UnsignedWord(int n = 0);
	UnsignedWord(const std::vector<long int> &n);
	UnsignedWord(const Bit &b);
	UnsignedWord(const UnsignedWord &u) : _bits(u._bits.size()) {
		for (int i = 0; i < bitLength(); ++i) {
			if (u._bits[i] != NULL)
				_bits[i] = new Bit(u[i]);
			else
				_bits[i] = NULL;
		}
	}

	~UnsignedWord() {
		for (int i = 0; i < bitLength(); ++i) {
			if (_bits[i] != NULL) {
				delete _bits[i];
				_bits[i] = NULL;
			}
		}
	}

	UnsignedWord<MAX_BIT_NUM, Bit> from_int(int c) {
		UnsignedWord<MAX_BIT_NUM, Bit> ret(c);
		return ret;
	}

	static UnsignedWord<MAX_BIT_NUM, Bit> static_from_int(int c) {
		UnsignedWord<MAX_BIT_NUM, Bit> ret(c);
		return ret;
	}

	static UnsignedWord<MAX_BIT_NUM, Bit> static_from_int(const std::vector<long int> &c) {
		UnsignedWord<MAX_BIT_NUM, Bit> ret(c);
		return ret;
	}

	void operator+=(const UnsignedWord &w);
	void operator-=(const UnsignedWord &w);
	void operator*=(const UnsignedWord &w) { *this = *this * w; }
	void operator*=(const Bit &b) { *this = *this * b; }

	void operator<<=(int i);
	void operator>>=(int i);

	Bit &operator[](int i) {
		assert(i < bitLength());
		if (_bits[i] == NULL)
			_bits[i] = new Bit;
		return *_bits[i];
	}
	const Bit &operator[](int i) const { return *_bits[i]; }

	//neg should expand the bit representation to MAX_BIT_NUM and use 2-complement
	//void neg();

	UnsignedWord<MAX_BIT_NUM, Bit> operator+(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c+=b; return c; }
	UnsignedWord<MAX_BIT_NUM, Bit> operator-(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c-=b; return c; }
	UnsignedWord<MAX_BIT_NUM, Bit> operator*(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const;
	UnsignedWord<MAX_BIT_NUM, Bit> operator*(const Bit &b) const;
	UnsignedWord<MAX_BIT_NUM, Bit> operator<<(int i) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c<<=i; return c; }
	UnsignedWord<MAX_BIT_NUM, Bit> operator>>(int i) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c>>=i; return c; }

	UnsignedWord<MAX_BIT_NUM, Bit> &operator=(const UnsignedWord<MAX_BIT_NUM, Bit> &u) {
		if (this == &u)
			return *this;

		for (unsigned int i = 0; i < _bits.size(); ++i) {
			if (_bits[i] != NULL) {
				delete _bits[i];
				_bits[i] = NULL;
			}
		}

		_bits.resize(u._bits.size());
		for (unsigned int i = 0; i < _bits.size(); ++i) {
			if (u._bits[i] != NULL) {
				_bits[i] = new Bit(*(u._bits[i]));
			}
		}
		return *this;
	}

	Bit operator==(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const;

	Bit operator>(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { return biggerBiggerEqual(b, false); }
	Bit operator<(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { return b.operator>(*this); }

	Bit operator>=(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { return biggerBiggerEqual(b, true); }
	Bit operator<=(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { return b.operator>=(*this); }

	Bit operator==(const int &b) const;
	Bit operator<(const int &b) const;
	Bit operator>(const int &b) const;

	int bitLength() const { return _bits.size(); }
	void set_bit_length(int len) {
		int oldLen = bitLength();

		if (len < bitLength()) {
			for (int i = len; i < bitLength(); ++i) {
				delete _bits[i];
			}
		}
		 _bits.resize(len);

		for (int i = oldLen; i < bitLength(); ++i)
			_bits[i] = NULL;
	}

	int size() const { return _bits.size(); }
	void set_bit(int i, const Bit &b) {
		assert(i < MAX_BIT_NUM);
		assert(i >= 0);
		if (_bits[i] != NULL)
			delete _bits[i];
		_bits[i] = new Bit(b);
	}
	unsigned long to_int() const;

	Bit bits_to_number() const {
		if (bitLength() == 0) {
			Bit ret(0);
			return ret;
		}
		Bit ret = (*this)[bitLength() - 1];

		for (int i = bitLength() - 2; i >= 0; --i) {
			ret += ret;
			ret += (*this)[i];
		}

		return ret;
	}

	std::string to_bit_stream() const {
		char buf[100];
		std::string ret = "";
		for (int i = bitLength() - 1; i >= 0; --i) {
			sprintf(buf, " %d", (*this)[i].to_int());
			ret += buf;
		}
		return ret;
	}

	template<int S, class B>
	friend std::ostream &operator<<(std::ostream &out, const UnsignedWord<S,B> &z);
	template<int S, class B>
	friend std::istream &operator>>(std::istream &in, UnsignedWord<S,B> &z);

	static int static_in_range(int a) { return a & (((unsigned int)-1) >> (32 - MAX_BIT_NUM)); }
	static int max_bit_num() { return MAX_BIT_NUM; }
};




template<int Size, class Bit>
inline std::ostream &operator<<(std::ostream &out, const UnsignedWord<Size, Bit> &z) {
	out << z._bits.size() << " ";

	for (unsigned int i = 0; i < z._bits.size(); ++i) {
		if (z._bits[i] != NULL) {
			out << "Y " << *(z._bits[i]);
		} else {
			out << "N ";
		}
	}

	return out;
}

template<int Size, class Bit>
inline std::istream &operator>>(std::istream &in, UnsignedWord<Size, Bit> &z) {
	int size;
	char c;

	in >> size;
	z.set_bit_length(size);

	for (unsigned int i = 0; i < z._bits.size(); ++i) {
		in >> c;
		if (c != 'N') {
			z._bits[i] = new Bit;
			in >> *(z._bits[i]);
		} else {
			if (z._bits[i] != NULL)
				delete z._bits[i];
			z._bits[i] = NULL;
		}
	}

	return in;
}


// implementation


template<int MAX_BIT_NUM, class Bit>
inline UnsignedWord<MAX_BIT_NUM, Bit>::UnsignedWord(int n) {
	while ((n > 0) && (_bits.size() < MAX_BIT_NUM)) {
		_bits.push_back(new Bit(n & 1));
		n >>= 1;
	}
	assert(n == 0);
}

template<int MAX_BIT_NUM, class Bit>
inline UnsignedWord<MAX_BIT_NUM, Bit>::UnsignedWord(const std::vector<long int>& _n) {

	std::vector<long int> n = _n;

	std::function<bool()> is_all_zero = [&n]()->bool {
		for (auto i = n.begin(); i != n.end(); ++i)
			if (*i != 0)
				return false;
		return true;
	};

	while ((_bits.size() < MAX_BIT_NUM) && (!is_all_zero())) {
		std::vector<long int> bit(n.size());
		for (unsigned int i = 0; i < n.size(); ++i)
			bit[i] = n[i] & 1;

		_bits.push_back(new Bit(bit));

		for (auto i = n.begin(); i != n.end(); ++i)
			(*i) >>= 1;
	}

	for (auto i = n.begin(); i != n.end(); ++i)
		assert(*i == 0);
}

template<int MAX_BIT_NUM, class Bit>
inline UnsignedWord<MAX_BIT_NUM, Bit>::UnsignedWord(const Bit &n) {
	_bits.push_back(new Bit(n));
}

template<int MAX_BIT_NUM, class Bit>
inline void UnsignedWord<MAX_BIT_NUM, Bit>::operator+=(const UnsignedWord<MAX_BIT_NUM, Bit> &w) {
	int len = std::max(bitLength(), w.bitLength());
	int min = std::min(bitLength(), w.bitLength());

	int i;
	Bit carry(0);

	// TODO: assumes p=2
	// specifically: a+b is assumed to be xor

	UnsignedWord<MAX_BIT_NUM, Bit> &bits = *this;

	for (i = 0; i < min; ++i) {
		Bit new_carry = w[i]*carry + bits[i]*(w[i] + carry);
		bits[i] += w[i] + carry;
		carry = new_carry;
	}

	for (i = min; i < bitLength(); ++i) {
		Bit new_carry =  bits[i]*carry;
		bits[i] += carry;
		carry = new_carry;
	}

	_bits.resize(std::min(len + 1, MAX_BIT_NUM));

	for (i = min; i < w.bitLength(); ++i) {
		Bit new_carry = w[i]*carry;
		bits[i] = w[i] + carry;
		carry = new_carry;
	}

	if (len < MAX_BIT_NUM)
		bits[len] = carry;
}

template<class Bit>
inline Bit Not(const Bit &b) {
	return Bit::static_from_int(1) - b;
}

template<int MAX_BIT_NUM, class Bit>
inline void UnsignedWord<MAX_BIT_NUM, Bit>::operator-=(const UnsignedWord<MAX_BIT_NUM, Bit> &w) {
	assert(w.bitLength() <= bitLength());

	int min = std::min(bitLength(), w.bitLength());

	int i;
	Bit borrow(0);
	UnsignedWord<MAX_BIT_NUM, Bit> &bits = *this;

	for (i = 0; i < min; ++i) {
		Bit new_borrow = Not(bits[i]) * (w[i] + borrow) + w[i]*borrow;
		bits[i] += w[i] + borrow;
		borrow = new_borrow;
	}

	for (i = min; i < bitLength(); ++i) {
		Bit new_borrow = Not(bits[i])*borrow;
		bits[i] += borrow;
		borrow = new_borrow;
	}
}

template<int MAX_BIT_NUM, class Bit>
void UnsignedWord<MAX_BIT_NUM, Bit>::operator>>=(int i) {
	if (i >= bitLength()) {
		set_bit_length(1);
		(*this)[0] = Bit(0);
		return;
	}
	if (i == 0)
		return;
	assert(i > 0);

	int j;
	for (j = 0; j < bitLength() - i; ++j) {
		if (_bits[j] != NULL)
			delete _bits[j];
		_bits[j] = _bits[j + i];
		_bits[j+i] = NULL;
	}
	set_bit_length(bitLength() - i);
}

template<int MAX_BIT_NUM, class Bit>
void UnsignedWord<MAX_BIT_NUM, Bit>::operator<<=(int i) {
	int j;
	if (i == 0)
		return;
	assert(i > 0);

	set_bit_length(std::min(MAX_BIT_NUM, bitLength() + i));

	for (j = bitLength() - 1; j >= i; --j) {
		if (_bits[j] != NULL)
			delete _bits[j];
		_bits[j] = _bits[j - i];
		_bits[j - i] = NULL;
	}
	for (j = i - 1; j >= 0; --j) {
		if (_bits[j] != NULL)
			delete _bits[j];
		_bits[j] = new Bit(0);
	}
}

template<int MAX_BIT_NUM, class Bit>
UnsignedWord<MAX_BIT_NUM, Bit> UnsignedWord<MAX_BIT_NUM, Bit>::operator*(const Bit &b) const {
	UnsignedWord<MAX_BIT_NUM, Bit> a(*this);
	for (int i = 0; i < bitLength(); ++i)
		a[i] *= b;
	return a;
}

template<int MAX_BIT_NUM, class Bit>
UnsignedWord<MAX_BIT_NUM, Bit> operator*(const Bit &b, const UnsignedWord<MAX_BIT_NUM, Bit> &w) { return w*b; }


template<int MAX_BIT_NUM, class Bit>
UnsignedWord<MAX_BIT_NUM, Bit> UnsignedWord<MAX_BIT_NUM, Bit>::operator*(const UnsignedWord<MAX_BIT_NUM, Bit> &w) const {
	if (bitLength() > w.bitLength())
		return w * (*this);

	BinomialTournament< UnsignedWord<MAX_BIT_NUM, Bit> > add(BinomialTournament< UnsignedWord<MAX_BIT_NUM, Bit> >::add);

	add.add_to_tournament((*this) * w[0]);

	for (int i = 1; i < w.bitLength(); ++i)
		add.add_to_tournament((*this << i) * w[i]);

	return add.unite_all();
}

template<int MAX_BIT_NUM, class Bit>
unsigned long UnsignedWord<MAX_BIT_NUM, Bit>::to_int() const {
	unsigned long r = 0;
	for (int i = 0; i < bitLength(); ++i)
		if (_bits[i] != NULL)
			r += (_bits[i]->to_int()) << i;
	return r;
}

// Assumes the bits of the nummber are wither 1 or 0
// output is a bit (i.e. 0 or 1)
// does not require the ring size to be 2
#define get_bit(x, w) (((x) >> (w)) & 1)
template<int MAX_BIT_NUM, class Bit>
Bit UnsignedWord<MAX_BIT_NUM, Bit>::operator<(const int &w) const {

	// nothing is smaller than 0 (in Z_p)
	if (w == 0)
		return Bit(0);

	// If our bitLength is shorter than bit length of w, we are smaller no need to continue
	if ((1 << bitLength()) - 1 < w)
		return Bit(1);

	// we know that bitLength() > wBitLength because previousl checked so

	BinomialTournament<Bit> isSmaller( BinomialTournament<Bit>::add );

	BinomialTournament<Bit> equalPrefix( BinomialTournament<Bit>::mul );
	for (int i = bitLength() - 1; i >= 0; --i) {
		if (get_bit(w, i) == 1) {
			isSmaller.add_to_tournament( equalPrefix.unite_all() * (Bit(1) - (*this)[i]) );
			equalPrefix.add_to_tournament( (*this)[i] );
		} else {
			equalPrefix.add_to_tournament( Bit(1) - (*this)[i] );
		}
	}

	return isSmaller.unite_all();
}


template<int MAX_BIT_NUM, class Bit>
Bit UnsignedWord<MAX_BIT_NUM, Bit>::operator>(const int &w) const {
	// If our bitLength is shorter than bit length of w, we are smaller no need to continue
	if ((1 << bitLength()) - 1 <= w) {
		return Bit(0);
	}

	// we know that bitLength() > wBitLength because previousl checked so

	BinomialTournament<Bit> isLarger( BinomialTournament<Bit>::add );

	BinomialTournament<Bit> equalPrefix( BinomialTournament<Bit>::mul );
	for (int i = bitLength() - 1; i >= 0; --i) {
		if (get_bit(w, i) == 0) {
			if (i == bitLength() - 1)
				isLarger.add_to_tournament( (*this)[i] );
			else
				isLarger.add_to_tournament( equalPrefix.unite_all() * (*this)[i] );
			equalPrefix.add_to_tournament( Bit(1) - (*this)[i] );
		} else {
			equalPrefix.add_to_tournament( (*this)[i] );
		}
	}

	return isLarger.unite_all();
}

template<int MAX_BIT_NUM, class Bit>
Bit UnsignedWord<MAX_BIT_NUM, Bit>::operator==(const int &w) const {
	// If we have the wrong number of bits we cant be equal
	if ((1 << bitLength()) - 1 < w) {
		return Bit(0);
	}

	// we know that bitLength() > wBitLength because previousl checked so

	BinomialTournament<Bit> equalPrefix( BinomialTournament<Bit>::mul );
	for (int i = bitLength() - 1; i >= 0; --i) {
		if (get_bit(w, i) == 0) {
			equalPrefix.add_to_tournament( Bit(1) - (*this)[i] );
		} else {
			equalPrefix.add_to_tournament( (*this)[i] );
		}
	}

	return equalPrefix.unite_all();
}
#undef get_bit

template<int MAX_BIT_NUM, class Bit>
Bit UnsignedWord<MAX_BIT_NUM, Bit>::biggerBiggerEqual(const UnsignedWord<MAX_BIT_NUM, Bit> &w, bool isBiggerEqual) const {
	if (bitLength() == 0)
		return Bit(0);
	if (w.bitLength() == 0)
		return Bit(1);

	BinomialTournament<Bit> samePrefix( BinomialTournament<Bit>::mul );
	BinomialTournament<Bit> biggerThan( BinomialTournament<Bit>::add );

	bool overZ2Field = ((*this)[0].get_ring_size() == 2);


	// take care of the prefix of *this (i.e. pad w with 0)
	if (bitLength() > w.bitLength()) {
		for (int bit = bitLength() - 1; bit > w.bitLength() - 1; --bit) {
			samePrefix.add_to_tournament( !((*this)[bit]) );
		}
		biggerThan.add_to_tournament( !(samePrefix.unite_all()) );
	}


	// take care of the prefix of w (i.e. pad *this with 0)
	if (w.bitLength() > bitLength()) {
		for (int bit = w.bitLength() - 1; bit > bitLength() - 1; --bit) {
			samePrefix.add_to_tournament( !w[bit] );
		}
	}

	// If the numbers are not the same bit-length treat is as if the smaller one is padded with zeros
	for (int bit = std::min(bitLength(), w.bitLength()) - 1; bit >= 0; --bit) {

		Bit m = (*this)[bit] * w[bit];

		if (samePrefix.is_empty()) {
			biggerThan.add_to_tournament( (*this)[bit] - m );
		} else {
			Bit samePref = samePrefix.unite_all();
			biggerThan.add_to_tournament( ((*this)[bit] - m) * samePref );
		}

		if ((bit > 0) || (isBiggerEqual)) {
			Bit sameBit = (*this)[bit] + w[bit];
			if (!overZ2Field) {
				sameBit -= m;
				sameBit -= m;
			}
			sameBit = !sameBit;
			samePrefix.add_to_tournament(sameBit);
		}
	}

	if (isBiggerEqual && (!samePrefix.is_empty()))
		biggerThan.add_to_tournament(samePrefix.unite_all());

	if (biggerThan.is_empty())
		return Bit(0);
	return biggerThan.unite_all();
}


template<int MAX_BIT_NUM, class Bit>
UnsignedWord<MAX_BIT_NUM, Bit> min(const UnsignedWord<MAX_BIT_NUM, Bit> &a, const UnsignedWord<MAX_BIT_NUM, Bit> &b) {
	Bit bIsGreater = a < b;
	return bIsGreater*a + (Bit(1)+bIsGreater)*b;
}

//template<int MAX_BIT_NUM, class Bit>
//std::ostream &operator<<(std::ostream &out, const UnsignedWord<Bit> &w) {
//	for (int i = w.bitLength() - 1; i >= 0; --i)
//		out << w[i];
//	return out;
//}

class TruncConversion {
public:
	template<int MAX_BIT_NUM, class Bit>
	static Bit convert(const UnsignedWord<MAX_BIT_NUM, Bit> &w) { return w[0]; }
};

#endif

