# Makefile for higher-order explicit demo
#
# Jacob Cadena
# Math 6321 @ SMU
# Fall 2016

CXX      = g++
CXXFLAGS = -O0 --std=c++11

HEADERS = matrix.hpp erk2.hpp erk3.hpp erk4.hpp ssp2.hpp ssp3.hpp ssp4.hpp weno5.hpp get_time.h

# executable targets
all : WENOConvergence.exe SSPConvergence.exe testStability_BR.exe testStability_BL.exe testStability_L.exe

WENOConvergence.exe	: WENOConvergence.o erk3.o weno5.o matrix.o get_time.o
	$(CXX) $(CXXFLAGS) -o $@ $^

SSPConvergence.exe : SSPConvergence.o ssp2.o ssp3.o ssp4.o matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^

testStability_BR.exe : testStability_BR.o ssp2.o ssp3.o ssp4.o erk2.o erk3.o erk4.o weno5.o matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^

testStability_BL.exe : testStability_BL.o ssp2.o ssp3.o ssp4.o erk2.o erk3.o erk4.o weno5.o matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^

testStability_L.exe : testStability_L.o ssp2.o ssp3.o ssp4.o erk2.o erk3.o erk4.o weno5.o matrix.o
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
