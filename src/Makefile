CCFLAGS = -g

NTLINCDIR = -I/home/hayim/academic/fhe/ntl-9.6.2/include
NTLLIBDIR = -L/home/hayim/academic/fhe/ntl-9.6.2/src

FHEINCDIR = -I/home/hayim/academic/fhe/HElib-master/src
FHELIBDIR = -L/home/hayim/academic/fhe/HElib-master/src

HEADSUPINCDIR = -I../include
INCS = $(NTLINCDIR) $(FHEINCDIR) $(HEADSUPINCDIR)

OBJ = helib_keys.o helib_number.o

libliphe.a: $(OBJ)
	ar rUv libliphe.a $(OBJ)

%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o libliphe.a
