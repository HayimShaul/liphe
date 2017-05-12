#ifndef ___COMPARE_LIPHE__
#define ___COMPARE_LIPHE__

template <class CLASS>
class CompareEuler {
private:
	CLASS &_val;
public:
	CompareEuler(CLASS &v) : _val(v) {}

	CLASS operator==(int b) const { return isZero_Euler(_val - b); }
	CLASS operator<=(int b) const;
	CLASS operator<(int b) const;
	CLASS operator>=(int b) const;
	CLASS operator>(int b) const;

	// similarly define other operators
};


#endif
