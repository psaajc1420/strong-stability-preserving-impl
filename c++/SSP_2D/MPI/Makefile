# Makefile for higher-order explicit demo
#
# Jacob Cadena
# Math 6321 @ SMU
# Fall 2016

CXX      = g++ -fopenmp
CXXFLAGS = -O2 --std=c++11
MPICXX = mpicxx

HEADERS = matrix.hpp weno_mpi.hpp erk3.hpp weno5.hpp get_time.h

# executable targets
all : WENOConvergence.exe 

WENOConvergence.exe : WENOConvergence.o weno_mpi.o erk3.o weno5.o matrix.o get_time.o
	$(MPICXX) $(CXXFLAGS) -o $@ $^

%.o : %.cpp $(HEADERS)
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

get_time.o : get_time.c get_time.h
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

# utilities
clean :
	\rm -rf *.o *.txt

realclean : clean
	\rm -rf *.exe *~
