#ifndef ___HEAAN_KEYS___
#define ___HEAAN_KEYS___

#include <Ring.h>
#include <SecretKey.h>
#include <Scheme.h>

class HeaanKeys {
private:
	Ring *_context;
	SecretKey *_secretKey;
	Scheme *_scheme;

	long _logQ;
	long _logP;
	long _logN;
public:
	HeaanKeys() : _context(NULL), _secretKey(NULL), _scheme(NULL), _logQ(0), _logP(0), _logN(0) {}
	HeaanKeys(const HeaanKeys &h) : _context(h._context), _secretKey(h._secretKey), _scheme(h._scheme), _logQ(h._logQ), _logP(h._logP), _logN(h._logN) {}

	void initKeys(long logN, long logQ, long logP) {
		_logN = logN;
		_logP = logP;
		_logQ = logQ;
		_context = new Ring(_logN, _logQ);
		_secretKey = new SecretKey(_context);
		_scheme = new Scheme(_secretKey, _context);
	}

	void encrypt(Ciphertext &c, float f) {
		complex<double> m;
		m.real(f);
		m.imag(0);

		c = _scheme->encrypt(&m, 1, _logP, _logQ);
	}

	double decrypt(Ciphertext &c) {
		complex<double> *m = _scheme->decrypt(_secretKey, &c);
		double ret = m->real();
		delete m;
		return ret;
	}

	Scheme *scheme() { return _scheme; }
	SecretKey *secretKey() { return _secretKey; }
	int logP() const { return _logP; }
	int logQ() const { return _logQ; }
	int logN() const { return _logN; }
};

#endif
