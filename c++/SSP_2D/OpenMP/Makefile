# Makefile for higher-order explicit demo
#
# Jacob Cadena
# Math 6321 @ SMU
# Fall 2016

CXX      = g++ -fopenmp
CXXFLAGS = -O2 --std=c++11

HEADERS = matrix.hpp erk3.hpp weno5.hpp get_time.h

# executable targets
all : WENOConvergence.exe 

WENOConvergence.exe : WENOConvergence.o erk3.o weno5.o matrix.o get_time.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

get_time.o : get_time.c get_time.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# utilities
clean :
	\rm -rf *.o *.txt

realclean : clean
	\rm -rf *.exe *~
