CCFLAGS = -g --std=c++11 -Wall

NTLINCDIR = -I../../ntl-10.5.0-multithread/include
NTLLIBDIR = -L../../ntl-10.5.0-multithread/src

HELIBINCDIR = -I../../HElib-multithread/src
HELIBLIBDIR = -L../../HElib-multithread/src

HEAANINCDIR = -I../../HEAAN_v2.1/HEAAN/src
HEAANLIBDIR = -L../../HEAAN_v2.1/HEAAN/src

LIPHEINCDIR = -I../include
INCS = $(NTLINCDIR) $(HELIBINCDIR) $(LIPHEINCDIR) $(HEAANINCDIR)

OBJ = zp.o primes.o thread_pool.o # heaan_number.o

all: libliphe.a libHelibLiphe.a

libliphe.a: $(OBJ)
	ar rUv $@ $^

libHelibLiphe.a: helib_keys.o helib_number.o
	ar rUv $@ $^

libHeaanLiphe.a: heaan_number.o
	ar rUv $@ $^

%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o *.a
