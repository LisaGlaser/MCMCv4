# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# I'm using CFlags for compiler options, and just assume that LDFLAGS works for libraries to include. I might be abusing it here
CFLAGS=-c
LDFLAGS=-I /media/glaser/data/Nottingham/Research/MCinNCG/Sims/SeriousProgramming/Eigen/Eigen3.24/

all: mcmc

mcmc: MCMCv5.o progParams.o Dirac.o
	$(CC) MCMCv5.o progParams.o Dirac.o -o MCMCv5

MCMCv5.o: MCMCv5.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) MCMCv5.cpp

progParams.o: progParams.cpp
	$(CC) $(CFLAGS) progParams.cpp

Dirac.o: Dirac.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) Dirac.cpp

clean:
	rm *.o MCMCv5
