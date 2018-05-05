/* Explicit 3th-order Runge-Kutta time stepper class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef ERK3_DEFINED__
#define ERK3_DEFINED__

// Inclusions
#include <math.h>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"


// Explicit RK3 time stepper class
class ERK3Stepper {

 private:

  RHSFunction *frhs;                        // pointer to ODE RHS function
  std::vector<double> z, f0, f1, f2;         // reused vectors
  Matrix A;                                 // Butcher table
  std::vector<double> b, c;
  
 public:

  // constructor (sets RHS function pointer, allocates local data)
  ERK3Stepper(RHSFunction& frhs_, std::vector<double>& y) { 
    frhs = &frhs_;                        // store RHSFunction pointer
    z  = y;   // allocate reusable data 
    f0 = y;   //   based on size of y
    f1 = y;
    f2 = y;
    A = Matrix(3,3);                      // Butcher table data
    A(1,0) =  0.5;
    A(2,0) = -1.0;
    A(2,1) =  2.0;
    b = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    c = {0.0, 0.5, 1.0};
  };

  // Evolve routine (evolves the solution)
  std::vector<double> Evolve(std::vector<double>& tspan, double h, std::vector<double>& y);

  // Single step calculation
  int Step(double t, double h, std::vector<double>& y);

};

#endif
