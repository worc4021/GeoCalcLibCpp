GMPROOT=/usr/local/Cellar/gmp/6.2.1
GMPCLASS=/Users/Manuel/Documents/Development/GeoCalcLibCpp/GmpCpp
EIGENROOT=/usr/local/Cellar/eigen/3.3.8_1/include/eigen3
LRSROOT=/Users/manuel/Documents/Development/MatlabPublic/GeoCalcLibCpp/lrslib
BOOSTROOT=/usr/local/Cellar/boost/1.74.0/include

OWNGMP=TRUE

IFLAGS=-I$(EIGENROOT) -I$(BOOSTROOT)/include -I$(LRSROOT) -I$(GMPROOT)/include

# CXX=g++ -v
# CC=gcc -v

LFLAGS=-L$(GMPROOT)/lib -lgmp
CFLAGS=-Wall

OBJ=lrsinterface.o lrslib.o lrsgmp.o

ifeq ($(OWNGMP),TRUE)
	IFLAGS+= -I$(GMPCLASS)
	LFLAGS+= -L$(GMPCLASS) -lgmppp
else
	CFLAGS+= -DGMPMATRIX
	OBJ+= gmpMatrix.o
	LFLAGS+= -lgmpxx
endif

CFLAGS+=$(IFLAGS)
LRSFLAGS=$(CFLAGS) -DGMP -DMP
CXXFLAGS=$(LRSFLAGS) -std=c++17


all: liblrsinterface.a

libgmppp.a: ../GmpCpp/mpz.cpp ../GmpCpp/mpq.cpp ../GmpCpp/mpqMatrix.cpp
	make -C ../GmpCpp -f Makefile libgmppp.a

liblrsinterface.a: $(OBJ)
	$(AR) -rcs liblrsinterface.a $(OBJ)

lrsgmp.o: $(LRSROOT)/lrsgmp.c
	$(CC) $(LRSFLAGS) -o lrsgmp.o -c $(LRSROOT)/lrsgmp.c 2> lrsgmp.err

lrslib.o: $(LRSROOT)/lrslib.c
	$(CC) $(LRSFLAGS) -o lrslib.o -c $(LRSROOT)/lrslib.c 2> lrslib.err

lrsinterface.o: lrsinterface.cpp gmpMatrix.hpp
	$(CXX) $(CXXFLAGS) -o lrsinterface.o -c lrsinterface.cpp 2> lrsinterface.err

gmpMatrix.o: gmpMatrix.cpp
	$(CXX) $(CXXFLAGS) -o gmpMatrix.o -c gmpMatrix.cpp 2> gmpMatrix.err

test.o: gmpcpptest.cpp
	$(CXX) $(CXXFLAGS) -o test.o -c gmpcpptest.cpp 2> test.err

test: test.o $(OBJ) libgmppp.a
	$(CXX) $(LFLAGS) -o test $(OBJ) test.o 2> test.log

clean:
	make -C ../GmpCpp -f Makefile clean
	rm $(OBJ) *.err *.a *.log