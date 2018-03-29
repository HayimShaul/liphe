CCFLAGS = -g --std=c++11 -Wall

NTLINCDIR = -I../../ntl-10.5.0-multithread/include
NTLLIBDIR = ../../ntl-10.5.0-multithread/src

FHEINCDIR = -I../../HElib-multithread/src
FHELIBDIR = ../../HElib-multithread/src

LIPHEINCDIR = -I../include
LIPHELIBDIR = -L../src

LIBS = $(LIPHELIBDIR) -lliphe $(FHELIBDIR)/fhe.a $(NTLLIBDIR)/ntl.a -lgmp -lpthread
INCS = $(NTLINCDIR) $(FHEINCDIR) $(LIPHEINCDIR)

all: test_zp test_helib

#all: test min fast_min min2 min3

test_zp: test_zp.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_helib: test_helib.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)


%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o test_zp test_helib

