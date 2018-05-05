/* Explicit 2nd-order Runge-Kutta time stepper class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef ERK2_DEFINED__
#define ERK2_DEFINED__

// Inclusions
#include <math.h>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"


// Explicit RK2 time stepper class
class ERK2Stepper {

 private:

  RHSFunction *frhs;   // pointer to ODE RHS function
  std::vector<double> f, z;    // reused vectors
  
 public:

  // constructor (sets RHS function pointer, allocates local data)
  ERK2Stepper(RHSFunction& frhs_, std::vector<double>& y) { 
    frhs = &frhs_;      // store RHSFunction pointer
    f  = y;             // allocate reusable data
    z  = y;             //   based on size of y
  };

  // Evolve routine (evolves the solution)
  std::vector<double> Evolve(std::vector<double>& tspan, double h, 
                             std::vector<double>& y);

};

#endif
