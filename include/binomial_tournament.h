#ifndef ___BINOMIAL_TOURNAMENT___
#define ___BINOMIAL_TOURNAMENT___

#include <utility>
#include <vector>

template<class Number>
class BinomialTournament {
private:
	typedef void UniteStrategy(Number &, const Number &);

	class IsEmptyBoolean {
	private:
		bool _is;
	public:
		IsEmptyBoolean() : _is(true) {}
		IsEmptyBoolean(bool b) : _is(b) {}
		IsEmptyBoolean(const BinomialTournament<Number>::IsEmptyBoolean &b) : _is(b._is) {}

		bool is() const { return _is; }
		IsEmptyBoolean operator=(bool b) { _is = b; return *this; }
		IsEmptyBoolean operator=(IsEmptyBoolean b) { _is = b._is; return *this; }
		operator bool() const { return is(); }
	};

	std::vector< std::pair<Number, IsEmptyBoolean> > _heap;

	void set_slot_empty(int i) { _heap[i].second = true; }
	void set_slot_full(int i) { _heap[i].second = false; }

	UniteStrategy *_unite_strategy;
public:
	BinomialTournament(UniteStrategy *unite_strategy) : _heap(10) { _unite_strategy = unite_strategy; }
	void add_to_tournament(const Number &, int level = 0);
	Number unite_all(UniteStrategy *u = NULL) const;
	int levels() const;
	int max_level() const { return _heap.size(); }
	bool is_slot_empty(int i) const { return _heap[i].second.is() ; }

	Number &number(int i) { return _heap[i].first; }
	const Number &number(int i) const { return _heap[i].first; }

	static void add(Number &x, const Number &y) { x += y; }
	static void mul(Number &x, const Number &y) { x *= y; }
};

template<class Number>
class AddBinomialTournament : public BinomialTournament<Number> {
public:
	AddBinomialTournament() : BinomialTournament<Number>(BinomialTournament<Number>::add) {}
};

template<class Number>
class MulBinomialTournament : public BinomialTournament<Number> {
public:
	MulBinomialTournament() : BinomialTournament<Number>(BinomialTournament<Number>::mul) {}
};

template<class Number>
int BinomialTournament<Number>::levels() const {
	int count = 0;
	for (int i = 0; i < _heap.size(); ++i)
		if (!is_slot_empty(i))
			++count;
	return count;
}

template<class Number>
void BinomialTournament<Number>::add_to_tournament(const Number &n, int level) {
	Number toAdd = n;

	while (1) {
		if (_heap.size() <= level) {
			_heap.resize(level + 1);
		}

		if (is_slot_empty(level)) {
			number(level) = toAdd;
			set_slot_full(level);
			return;
		} else {
			_unite_strategy(toAdd, number(level));
			set_slot_empty(level);
			++level;
		}
	}
}

template<class Number>
Number BinomialTournament<Number>::unite_all(UniteStrategy *unite_strategy) const {
	if (unite_strategy == NULL)
		unite_strategy = _unite_strategy;

	Number n;
	int i = 0;
	while ((i < _heap.size()) && is_slot_empty(i))
		++i;

	if (i < _heap.size())
		n = number(i);
	++i;

	while (i < _heap.size()) {
		if (!is_slot_empty(i)) {
			unite_strategy(n, number(i));
		}
		++i;
	}

	return n;
}

#endif
