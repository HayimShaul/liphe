CCFLAGS = -g --std=c++11 -Wall

NTLINCDIR = -I/home/hayim/academic/fhe/ntl-9.6.2/include
NTLLIBDIR = -L/home/hayim/academic/fhe/ntl-9.6.2/src

FHEINCDIR = -I/home/hayim/academic/fhe/HElib-master/src
FHELIBDIR = /home/hayim/academic/fhe/HElib-master/src

LIPHEINCDIR = -I../include
LIPHELIBDIR = -L../src

LIBS = $(LIPHELIBDIR) -lliphe $(FHELIBDIR)/fhe.a $(NTLLIBDIR) -lntl
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

