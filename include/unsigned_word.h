#ifndef ___UNSIGNED_WORD___
#define ___UNSIGNED_WORD___

#include <assert.h>
#include <bitset>
#include <vector>

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
	std::vector<Bit> _bits;
public:
	UnsignedWord(int n = 0);
	UnsignedWord(const Bit &b);

	~UnsignedWord() {
	}

	UnsignedWord<MAX_BIT_NUM, Bit> from_int(int c) {
		UnsignedWord<MAX_BIT_NUM, Bit> ret(c);
		return ret;
	}

	static UnsignedWord<MAX_BIT_NUM, Bit> static_from_int(int c) {
		UnsignedWord<MAX_BIT_NUM, Bit> ret(c);
		return ret;
	}

	void operator+=(const UnsignedWord &w);
	void operator-=(const UnsignedWord &w);
	void operator*=(const UnsignedWord &w) { *this = *this * w; }
	void operator*=(const Bit &b) { *this = *this * b; }

	void operator<<=(int i);
	void operator>>=(int i);

	Bit &operator[](int i) { return _bits[i]; }
	const Bit &operator[](int i) const { return _bits[i]; }

	//neg should expand the bit representation to MAX_BIT_NUM and use 2-complement
	//void neg();

	UnsignedWord<MAX_BIT_NUM, Bit> operator+(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c+=b; return c; }
	UnsignedWord<MAX_BIT_NUM, Bit> operator-(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c-=b; return c; }
	UnsignedWord<MAX_BIT_NUM, Bit> operator*(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const;
	UnsignedWord<MAX_BIT_NUM, Bit> operator*(const Bit &b) const;
	UnsignedWord<MAX_BIT_NUM, Bit> operator<<(int i) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c<<=i; return c; }
	UnsignedWord<MAX_BIT_NUM, Bit> operator>>(int i) const { UnsignedWord<MAX_BIT_NUM, Bit> c(*this); c>>=i; return c; }

	Bit operator>(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const;
	Bit operator<(const UnsignedWord<MAX_BIT_NUM, Bit> &b) const { return b.operator>(*this); }

	Bit operator<(const int &b) const;
	Bit operator>(const int &b) const;

//	void setBitLength(int i) { _bits.resize(i); }
	int bitLength() const { return _bits.size(); }
	int size() const { return _bits.size(); }
	unsigned long to_int() const;
//	std::ostream &operator<<(std::ostream &out);

	static int static_in_range(int a) { return a & (((unsigned int)-1) >> (32 - MAX_BIT_NUM)); }
	static int max_bit_num() { return MAX_BIT_NUM; }
};

// implementation


template<int MAX_BIT_NUM, class Bit>
inline UnsignedWord<MAX_BIT_NUM, Bit>::UnsignedWord(int n) {
	while ((n > 0) && (_bits.size() < MAX_BIT_NUM)) {
		_bits.push_back(Bit(n & 1));
		n >>= 1;
	}
}

template<int MAX_BIT_NUM, class Bit>
inline UnsignedWord<MAX_BIT_NUM, Bit>::UnsignedWord(const Bit &n) {
	_bits.push_back(n);
}

template<int MAX_BIT_NUM, class Bit>
inline void UnsignedWord<MAX_BIT_NUM, Bit>::operator+=(const UnsignedWord<MAX_BIT_NUM, Bit> &w) {
	int len = std::max(bitLength(), w.bitLength());
	int min = std::min(bitLength(), w.bitLength());

	int i;
	Bit carry(0);

	for (i = 0; i < min; ++i) {
		Bit new_carry = w[i]*carry + _bits[i]*(w[i] + carry);
		_bits[i] += w[i] + carry;
		carry = new_carry;
	}

	for (i = min; i < bitLength(); ++i) {
		Bit new_carry =  _bits[i]*carry;
		_bits[i] += carry;
		carry = new_carry;
	}

	_bits.resize(std::min(len + 1, MAX_BIT_NUM));

	for (i = min; i < w.bitLength(); ++i) {
		Bit new_carry = w[i]*carry;
		_bits[i] = w[i] + carry;
		carry = new_carry;
	}

	if (len < MAX_BIT_NUM)
		_bits[len] = carry;
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

	for (i = 0; i < min; ++i) {
		Bit new_borrow = Not(_bits[i]) * (w[i] + borrow) + w[i]*borrow;
		_bits[i] += w[i] + borrow;
		borrow = new_borrow;
	}

	for (i = min; i < bitLength(); ++i) {
		Bit new_borrow = Not(_bits[i])*borrow;
		_bits[i] += borrow;
		borrow = new_borrow;
	}
}

template<int MAX_BIT_NUM, class Bit>
void UnsignedWord<MAX_BIT_NUM, Bit>::operator>>=(int i) {
	if (i >= bitLength()) {
		_bits.resize(1);
		_bits[0] = Bit(0);
		return;
	}

	int j;
	for (j = 0; j < bitLength() - i; ++j)
		_bits[j] = _bits[j + i];
	_bits.resize(bitLength() - i);
}

template<int MAX_BIT_NUM, class Bit>
void UnsignedWord<MAX_BIT_NUM, Bit>::operator<<=(int i) {
	int j;

	_bits.resize(std::min(MAX_BIT_NUM, bitLength() + i));

	for (j = bitLength() - 1; j >= i; --j)
		_bits[j] = _bits[j - i];
	for (j = i - 1; j >= 0; --j)
		_bits[j] = Bit(0);
}

template<int MAX_BIT_NUM, class Bit>
UnsignedWord<MAX_BIT_NUM, Bit> UnsignedWord<MAX_BIT_NUM, Bit>::operator*(const Bit &b) const {
	UnsignedWord<MAX_BIT_NUM, Bit> a(*this);
	for (int i = 0; i < bitLength(); ++i)
		a._bits[i] *= b;
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
		r += (_bits[i].to_int()) << i;
	return r;
}

// Assumes the bits of the nummber are wither 1 or 0
// output is a bit (i.e. 0 or 1)
// does not require the ring size to be 2
#define get_bit(x, w) (((x) >> (w)) & 1)
template<int MAX_BIT_NUM, class Bit>
Bit UnsignedWord<MAX_BIT_NUM, Bit>::operator<(const int &w) const {
std::cerr << "UnignedWord:::operator<" << std::endl;

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
	if ((1 << bitLength()) - 1 < w) {
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
			equalPrefix.add_to_tournament( (*this)[i] );
		} else {
			equalPrefix.add_to_tournament( Bit(1) - (*this)[i] );
		}
	}

	return isLarger.unite_all();
}
#undef get_bit

template<int MAX_BIT_NUM, class Bit>
Bit UnsignedWord<MAX_BIT_NUM, Bit>::operator>(const UnsignedWord<MAX_BIT_NUM, Bit> &w) const {
	int i;

	if (bitLength() == 0)
		return Bit(0);
	if (w.bitLength() == 0)
		return Bit(1);

	Bit ourMsbIsAllZero(1);
	Bit hisMsbIsAllZero(1);

	//	compute this loop in a shallower circuit
	//	for (i = bitLength() - 1; i >= w.bitLength(); --i)
	//		ourMsbIsAllZero *= !((*this)[i]);
	if (bitLength() > w.bitLength()) {
		std::vector<Bit> arr;
		arr.resize(bitLength() - w.bitLength());
		for (int j = w.bitLength(); j < bitLength(); ++j)
			arr[j - w.bitLength()] = !((*this)[j]);
		ourMsbIsAllZero = mulArray(arr);
	}

	//	compute this loop in a shallower circuit
	//	for (i = w.bitLength() - 1; i >= bitLength(); --i)
	//		hisMsbIsAllZero *= !(w[i]);
	if (w.bitLength() > bitLength()) {
		std::vector<Bit> arr;
		arr.resize(w.bitLength() - bitLength());
		for (int j = bitLength(); j < w.bitLength(); ++j)
			arr[j - bitLength()] = !(w[j]);
		hisMsbIsAllZero = mulArray(arr);
	}

	int commonLength = std::min(bitLength(), w.bitLength());

	std::vector<Bit> sameBit;
	sameBit.resize(commonLength);
	for (i = 0; i < commonLength; ++i)
		sameBit[i] = !((*this)[i] + w[i]);

	std::vector<Bit> equalMsbVec;
	equalMsbVec.resize(commonLength);
	equalMsbVec[commonLength - 1] = Bit(1);
	for (i = 0; i < commonLength - 1; ++i)
		equalMsbVec[i] = mulArray(sameBit, i+1, -1);

// TODO: bodyIsGreaterVec should be bitDeterminesGreater or something
	std::vector<Bit> bodyIsGreaterVec;
	bodyIsGreaterVec.resize(commonLength);
	for (i = 0; i < commonLength - 1; ++i)
		bodyIsGreaterVec[i] = equalMsbVec[i]*((*this)[i])*(!w[i]);
	// We know equalMsbVec[commonLength-1] == Bit(1)
	bodyIsGreaterVec[commonLength - 1] = ((*this)[commonLength - 1])*(!w[commonLength - 1]);
	Bit bodyIsGreater = addArray(bodyIsGreaterVec);

//	Bit equalMSB(1);
//	Bit bodyIsGreater(0);
//	for (i = commonLength - 1; i >= 0; --i) {
////		bodyIsGreater = bodyIsGreater | (equalMSB*((*this)[i])*(!w[i]));
////		which equals
////		bodyIsGreater = bodyIsGreater + (equalMSB*((*this)[i])*(!w[i])) + bodyIsGreater*(equalMSB*((*this)[i])*(!w[i]));
////		but bodyIsGreater==1 =>  equalMSB==0,    and equalMSB==1 => bodyIsGreater==0,   so we can write:
////		bodyIsGreater = bodyIsGreater + (equalMSB*((*this)[i])*(!w[i]));
////		equalMSB *= !(((*this)[i]) + w[i]);
//		bodyIsGreater = bodyIsGreater + (equalMsbVec[i]*((*this)[i])*(!w[i]));
//	}

//	Bit isGreater = (!ourMsbIsAllZero) | (hisMsbIsAllZero*bodyIsGreater);
//	Bit isGreater = (!ourMsbIsAllZero) + (hisMsbIsAllZero*bodyIsGreater) + (!ourMsbIsAllZero)*hisMsbIsAllZero*bodyIsGreater;
//	Since ourMsbAllZero==0  => hisMsbIsAllZero==1  we can write
//	Bit isGreater = (!ourMsbIsAllZero) + (hisMsbIsAllZero*bodyIsGreater) + (!ourMsbIsAllZero)*bodyIsGreater;
//	but that is really:
//	Bit isGreater = (Bit(1) + ourMsbIsAllZero) + (hisMsbIsAllZero*bodyIsGreater) + (Bit(1) + ourMsbIsAllZero)*bodyIsGreater;
//	Bit isGreater = Bit(1) + ourMsbIsAllZero + (hisMsbIsAllZero*bodyIsGreater) + (ourMsbIsAllZero*bodyIsGreater) + bodyIsGreater;
//	Bit isGreater = Bit(1) + ourMsbIsAllZero + (Bit(1)+hisMsbIsAllZero+ourMsbIsAllZero)*bodyIsGreater;
	Bit isGreater = Bit(1) + ourMsbIsAllZero + (Bit(1)+hisMsbIsAllZero+ourMsbIsAllZero)*bodyIsGreater;
//	Bit isGreater = Bit(1) + ourMsbIsAllZero + (!(hisMsbIsAllZero+ourMsbIsAllZero))*bodyIsGreater;

	return isGreater;
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

