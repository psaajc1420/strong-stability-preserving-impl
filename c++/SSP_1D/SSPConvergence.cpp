/* Main routine to test convergence of 2nd, 3rd and 4th order
   Strong Stability Preserving Methods.

     y' = -y, t in [0,1],
     y(0) = 1.

   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "ssp2.hpp"
#include "ssp3.hpp"
#include "ssp4.hpp"
#include "weno5.hpp"


using namespace std;


// ODE RHS function
class MyRHS: public RHSFunction {
public:
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = -y[0];
    return 0;
  }
};

// Convenience function for analytical solution
vector<double> ytrue(const double t) { 
  vector<double> yt = {exp(-t)};
  return yt;
};


// main routine
int main() {

  // time steps to try
  vector<double> Nt = {20, 40, 80, 160, 320};

  // set problem information
  vector<double> y0 = {1.0};
  double t0 = 0.0;
  double Tf = 1.0;

  // create ODE RHS function objects
  MyRHS rhs;

  // create SSP stepper objects
  SSP2Stepper SSP2(rhs, y0);
  SSP3Stepper SSP3(rhs, y0);
  SSP4Stepper SSP4(rhs, y0);

  // storage for error results
  vector<double> errs(Nt.size());

  ///////// SSP2 /////////
  cout << "\nSSP2 Method:\n";

  // loop over time step sizes
  for (int ih=0; ih<Nt.size(); ih++) {

    double dtout = (Tf-t0)/Nt[ih];
    double h = dtout;
    // set the initial condition, initial time
    vector<double> y(y0);
    double tcur = t0;

    // reset maxerr
    double maxerr = 0.0;
    
    // loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt[ih]; istep++) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};
 
      // call the solver, update current time
      vector<double> tvals = SSP2.Evolve(tspan, h, y);
      tcur = tvals.back();     // last entry in tvals

      // calculate exact value at each time
      vector<double> exact = ytrue(tcur);

      // calculate the error
      vector<double> yerr = y - exact;
      double err = InfNorm(yerr);
      maxerr = std::max(maxerr, err);

    }
    
    // compute convergence rate
    cout << "   n = " << Nt[ih] << "\t  error = " << maxerr;
    errs[ih] = maxerr;
    if (ih > 0)
      cout << "\t  conv rate = " << -(log(errs[ih])-log(errs[ih-1]))/(log(Nt[ih])-log(Nt[ih-1]));
    cout << endl;


  }
    Matrix error(1, Nt.size());
    Matrix n(1,Nt.size());
    for(int j = 0; j < Nt.size(); j++)
      n[j] = {Nt[j]};
    for(int j = 0; j < Nt.size(); j++)
      error[j] = {errs[j]};
    n.Write("ntvals_ssp2.txt");
    error.Write("error_ssp2.txt");


///////// SSP3 /////////
  cout << "\nSSP3 Method:\n";

  // loop over time step sizes
  for (int ih=0; ih<Nt.size(); ih++) {

    double dtout = (Tf-t0)/Nt[ih];
    double h = dtout;
    // set the initial condition, initial time
    vector<double> y(y0);
    double tcur = t0;

    // reset maxerr
    double maxerr = 0.0;
    
    // loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt[ih]; istep++) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};
 
      // call the solver, update current time
      vector<double> tvals = SSP3.Evolve(tspan, h, y);
      tcur = tvals.back();     // last entry in tvals

      // calculate exact value at each time
      vector<double> exact = ytrue(tcur);

      // calculate the error
      vector<double> yerr = y - exact;
      double err = InfNorm(yerr);
      maxerr = std::max(maxerr, err);

    }
    
    // compute convergence rate
    cout << "   n = " << Nt[ih] << "\t  error = " << maxerr;
    errs[ih] = maxerr;
    if (ih > 0)
      cout << "\t  conv rate = " << -(log(errs[ih])-log(errs[ih-1]))/(log(Nt[ih])-log(Nt[ih-1]));
    cout << endl;


  }
    for(int j = 0; j < Nt.size(); j++)
      error[j] = {errs[j]};
    n.Write("ntvals_ssp3.txt");
    error.Write("error_ssp3.txt");




 ///////// SSP4 /////////
  cout << "\nSSP4 Method:\n";

  // loop over time step sizes
  for (int ih=0; ih<Nt.size(); ih++) {

    double dtout = (Tf-t0)/Nt[ih];
    double h = dtout;
    // set the initial condition, initial time
    vector<double> y(y0);
    double tcur = t0;

    // reset maxerr
    double maxerr = 0.0;
    
    // loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt[ih]; istep++) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};
 
      // call the solver, update current time
      vector<double> tvals = SSP4.Evolve(tspan, h, y);
      tcur = tvals.back();     // last entry in tvals

      // calculate exact value at each time
      vector<double> exact = ytrue(tcur);

      // calculate the error
      vector<double> yerr = y - exact;
      double err = InfNorm(yerr);
      maxerr = std::max(maxerr, err);

    }
    
    // compute convergence rate
    cout << "   n = " << Nt[ih] << "\t  error = " << maxerr;
    errs[ih] = maxerr;
    if (ih > 0)
      cout << "\t  conv rate = " << -(log(errs[ih])-log(errs[ih-1]))/(log(Nt[ih])-log(Nt[ih-1]));
    cout << endl;


  }
    for(int j = 0; j < Nt.size(); j++)
      error[j] = {errs[j]};
    n.Write("ntvals_ssp4.txt");
    error.Write("error_ssp4.txt");


  return 0;
}
