#ifndef ___COMPARE_LIPHE__
#define ___COMPARE_LIPHE__

template <class Number>
class CompareEuler {
private:
	const Number &_val;
public:
	CompareEuler(const Number &v) : _val(v) {}

	Number operator==(int b) const { return isZero_Euler(((b == 0) ? _val : (_val - b))); }
	Number operator!=(int b) const { return Number(1) - operator==(b); }
	Number operator<=(int b) const;
	Number operator<(int b) const;
	Number operator>=(int b) const;
	Number operator>(int b) const;
};

template <class Number>
class CompareNative {
private:
	const Number &_val;
public:
	CompareNative(const Number &v) : _val(v) {}

	Number operator==(int b) const { return _val == b; }
	Number operator<=(int b) const { return _val <= b; }
	Number operator<(int b) const { return _val < b; }
	Number operator>=(int b) const { return _val >= b; }
	Number operator>(int b) const { return _val > b; }
};

#endif
