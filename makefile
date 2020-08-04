# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# I'm using CFlags for compiler options, and just assume that LDFLAGS works for libraries to include. I might be abusing it here
CFLAGS= -c -IEigen/Eigen3.24 -g3 -std=c++11 -Wno-ignored-attributes -Wno-deprecated-declarations -IEigen/Eigen3.24 
# LDFLAGS= -I /usr/local/include/eigen3
LDFLAGS= -I eigen-3.2.4

all: mcmc

mcmc: MCMCv4.o progParams.o Dirac.o
	$(CC) MCMCv4.o progParams.o Dirac.o -o MCMCv4

MCMCv4.o: MCMCv4.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) MCMCv4.cpp

progParams.o: progParams.cpp
	$(CC) $(CFLAGS) progParams.cpp

Dirac.o: Dirac.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) Dirac.cpp

clean:
	rm *.o MCMCv4
