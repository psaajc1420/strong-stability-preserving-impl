/* Explicit 3th-order Runge-Kutta time stepper class implementation file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#include "erk3.hpp"



std::vector<double> ERK3Stepper::Evolve(std::vector<double>& tspan, double h, 
                                        std::vector<std::vector<double>> &y) {

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

  int i;

  // iterate over time steps
  for (i=0; i<N; i++) {

    // last step only: update h to stop directly at final time
    if (i == N-1) 
      h = tspan[1]-times[i];

    // perform a single step of RK3 to update y
    if (Step(times[i], h, y) != 0) {
      std::cerr << "Evolve: Error in Step() function\n";
    }
    // update current time, store in output array
    times.push_back(times[i] + h);
  }

  return times;
}


// Single step of explicit 3th-order Runge-Kutta
//
// Inputs:  t holds the current time
//          h holds the current time step size
//          z, f1-f3 hold temporary vectors needed for the problem
//          y holds the current solution
// Outputs: y holds the updated solution, y(t+h)
//
// The return value is an integer indicating success/failure,
// with 0 indicating success, and nonzero failure.
int ERK3Stepper::Step(double t, double h, std::vector<std::vector<double>> &y) {
 
  // Intialize local variables
  int i, j;
  
  // stage 1: set stage and compute RHS
  z = y;
  if (frhs->Evaluate(t, z, f0) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }
  
  // stage 2: set stage and compute RHS
  // # pragma omp parallel for collapse(2) private(i,j)
  for(i = 0; i < y.size(); i++)
    for(j = 0; j < y[0].size(); j++)
      z[i][j] = y[i][j] + h*A(1,0)*f0[i][j];

  if (frhs->Evaluate(t+c[1]*h, z, f1) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }
   
  // stage 3: set stage and compute RHS
  // # pragma omp parallel for collapse(2) private(i,j)
  for(i = 0; i < y.size(); i++)
    for(j = 0; j < y[0].size(); j++)
      z[i][j] = y[i][j] + h*(A(2,0)*f0[i][j] + A(2,1)*f1[i][j]);

  if (frhs->Evaluate(t+c[2]*h, z, f2) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // compute next step solution
  // # pragma omp parallel for collapse(2) private(i,j)
  for(i = 0; i < y.size(); i++)
    for(j = 0; j < y[0].size(); j++)
      y[i][j] = y[i][j] +  h*(b[0]*f0[i][j] + b[1]*f1[i][j] + b[2]*f2[i][j]);

  // return success
  return 0;
}
