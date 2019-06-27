CCFLAGS = -g --std=c++11 -Wall

NTLINCDIR = -I../../ntl-10.5.0-multithread/include
NTLLIBDIR = -L../../ntl-10.5.0-multithread/src

HELIBINCDIR = -I../../HElib-multithread/src
HELIBLIBDIR = -L../../HElib-multithread/src

HEAANINCDIR = -I../../HEAAN/HEAAN/src
HEAANLIBDIR = -L../../HEAAN/HEAAN/src

LIPHEINCDIR = -I../include
INCS = $(NTLINCDIR) $(HELIBINCDIR) $(LIPHEINCDIR) $(HEAANINCDIR)

OBJ = primes.o thread_pool.o # heaan_number.o

all: libliphe.a libHelibLiphe.a

libliphe.a: $(OBJ)
	ar rUv libliphe.a $(OBJ)

libHelibLiphe.a: helib_keys.o helib_number.o
	ar rUv libHelibLiphe.a $(OBJ)

libHeaanLiphe.a: heaan_number.o
	ar rUv libHeaanLiphe.a $(OBJ)

%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o libliphe.a
