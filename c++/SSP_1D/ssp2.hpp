/* Explicit 2nd-order SSP time stepper class header file.

   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

#ifndef SSP2_DEFINED__
#define SSP2_DEFINED__

// Inclusions
#include <math.h>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"


// Explicit SSP2 time stepper class
class SSP2Stepper {

 private:

  RHSFunction *frhs;                      // pointer to ODE RHS function
  std::vector<double> z, f0, f1;          // reused vectors
  
 public:

  // constructor (sets RHS function pointer, allocates local data)
  SSP2Stepper(RHSFunction& frhs_, std::vector<double>& y) { 
    frhs = &frhs_;                        // store RHSFunction pointer
    z  = y;   // allocate reusable data 
    f0 = y;   //   based on size of y
    f1 = y;
  };

  // Evolve routine (evolves the solution)
  std::vector<double> Evolve(std::vector<double>& tspan, double h, std::vector<double>& y);

  // Single step calculation
  int Step(double t, double h, std::vector<double>& y);

};

#endif
