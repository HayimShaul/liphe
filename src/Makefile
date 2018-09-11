CCFLAGS = -g --std=c++11 -Wall

NTLINCDIR = -I../../ntl-10.5.0/include
NTLLIBDIR = -L../../ntl-10.5.0/src

HELIBINCDIR = -I../../HElib/src
HELIBLIBDIR = -L../../HElib/src

HEAANINCDIR = -I../../HEAAN/HEAAN/src
HEAANLIBDIR = -L../../HEAAN/HEAAN/src

LIPHEINCDIR = -I../include
INCS = $(NTLINCDIR) $(HELIBINCDIR) $(LIPHEINCDIR) $(HEAANINCDIR)

OBJ = helib_keys.o helib_number.o primes.o thread_pool.o heaan_number.o

libliphe.a: $(OBJ)
	ar rUv libliphe.a $(OBJ)

%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o libliphe.a
