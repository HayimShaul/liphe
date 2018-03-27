CCFLAGS = -g --std=c++11 -Wall

NTLINCDIR = -I../../ntl-10.5.0/include
NTLLIBDIR = -L../../ntl-10.5.0/src

FHEINCDIR = -I../../HElib/src
FHELIBDIR = -L../../HElib/src

HEADSUPINCDIR = -I../include
INCS = $(NTLINCDIR) $(FHEINCDIR) $(HEADSUPINCDIR)

OBJ = helib_keys.o helib_number.o primes.o thread_pool.o

libliphe.a: $(OBJ)
	ar rUv libliphe.a $(OBJ)

%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o libliphe.a
