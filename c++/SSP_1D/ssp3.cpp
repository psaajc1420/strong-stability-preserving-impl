/* Explicit 3rd-order SSP time stepper class implementation file.
  
   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

#include <vector>
#include "matrix.hpp"
#include "ssp3.hpp"


// The explicit 3rd-order SSP time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          h holds the desired time step size
//          y holds the initial condition, y(t0)
// Outputs: y holds the computed solution, y(tf)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> SSP3Stepper::Evolve(std::vector<double>& tspan, double h, 
                                        std::vector<double>& y) {

  // initialize output
  std::vector<double> times = {tspan[0]};

  // check for legal inputs 
  if (h <= 0.0) {
    std::cerr << "Evolve: Illegal h\n";
    return times;
  }
  if (tspan[1] <= tspan[0]) {
    std::cerr << "Evolve: Illegal tspan\n";
    return times;	  
  }
  

  // figure out how many time steps
  long int N = (tspan[1]-tspan[0])/h;
  if (tspan[1] > tspan[0]+N*h)  N++;

  // iterate over time steps
  for (int i=0; i<N; i++) {

    // last step only: update h to stop directly at final time
    if (i == N-1) 
      h = tspan[1]-times[i];

    // perform a single step of RK3 to update y
    if (Step(times[i], h, y) != 0) {
      std::cerr << "Evolve: Error in Step() function\n";
      return times;
    }
    // update current time, store in output array
    times.push_back(times[i] + h);
  }

  return times;
}


// Single step of explicit 3rd-order SSP
//
// Inputs:  t holds the current time
//          h holds the current time step size
//          z, f1-f3 hold temporary vectors needed for the problem
//          y holds the current solution
// Outputs: y holds the updated solution, y(t+h)
//
// The return value is an integer indicating success/failure,
// with 0 indicating success, and nonzero failure.
int SSP3Stepper::Step(double t, double h, std::vector<double>& y) {

  // stage 1: set stage and compute RHS
  z = y;
  if (frhs->Evaluate(t, z, f0) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 2: set stage and compute RHS
  z = y + h*f0;
  if (frhs->Evaluate(t, z, f1) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 3: set stage and compute RHS
  z = 3.0/4.0*y + 1.0/4.0*(z + h*f1);
  if (frhs->Evaluate(t, z, f2) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  } 
  //compute next step solution
  y = (y + 2.0*z + 2.0*h*f2)/3.0;

  // return success
  return 0;
}
