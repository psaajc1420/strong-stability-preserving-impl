/* Explicit 4th-order Runge-Kutta time stepper class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef ERK4_DEFINED__
#define ERK4_DEFINED__

// Inclusions
#include <math.h>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"


// Explicit RK4 time stepper class
class ERK4Stepper {

 private:

  RHSFunction *frhs;                        // pointer to ODE RHS function
  std::vector<double> z, f0, f1, f2, f3;    // reused vectors
  Matrix A;                                 // Butcher table
  std::vector<double> b, c;

 public:

  // constructor (sets RHS function pointer, allocates local data)
  ERK4Stepper(RHSFunction& frhs_, std::vector<double>& y) { 
    frhs = &frhs_;                        // store RHSFunction pointer
    z  = y;   // allocate reusable data 
    f0 = y;   //   based on size of y
    f1 = y;
    f2 = y;
    f3 = y;
    A = Matrix(4,4);                      // Butcher table data
    A(1,0) = 0.5;
    A(2,1) = 0.5;
    A(3,2) = 1.0;
    b = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    c = {0.0, 0.5, 0.5, 1.0};
  };

  // Evolve routine (evolves the solution)
  std::vector<double> Evolve(std::vector<double>& tspan, double h, std::vector<double>& y);

  // Single step calculation
  int Step(double t, double h, std::vector<double>& y);

};

#endif
