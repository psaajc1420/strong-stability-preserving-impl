/* Main routine to test convergence of fifth order WENO method with 
   the a linear flux

     y' = L(u), t in [0,1],     // ODE
     L(u) = u                   // Burgers Flux
     y(0) = sin(pi*x)           // Initial Condition
   
   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "erk3.hpp"
#include "weno5.hpp"
#include "get_time.h"


using namespace std;


// ODE RHS function
class MyRHS: public RHSFunction {
private:
  WEN05 *weno;
public:
  
  MyRHS(std::vector<double>& u_, RHSFlux& flux_, RHSSource& source_, double dx_){
    weno = new WEN05(u_, flux_, source_, dx_);
  }
  int Evaluate(double t, vector<double>& y, vector<double>& f) { 
    f = weno->residual(y);
    return 0;
  }
  ~MyRHS(){
    delete weno;
  }
};

// ODE Flux function
class MyFlux: public RHSFlux {
public:
  int Evaluate(vector<double>& y, vector<double>& f, vector<double>& df) {
    for(int i = 0; i < y.size(); i++)
      f[i] = y[i];

    for(int i = 0; i < y.size(); i++)
      df[i] = 1;

    return 0;
  }
};

// ODE Source function
class MySource : public RHSSource{
public:
  int Evaluate(vector<double>& y, vector<double>& f) {
    for(int i = 0; i  < y.size(); i++)
      f[i] = 0;
    return 0;
  }
};

// Convenience function for analytical solution
vector<double> ytrue(const double t, vector<double> x) { 
  vector<double> yt;
  for(int i = 0; i  < x.size(); i++)
    yt.push_back(sin(M_PI*(x[i] - t)));
  return yt;
};


// main routine
int main() {

  // time steps to try
<<<<<<< HEAD
  // vector<double> nx = {1200};
=======
>>>>>>> master
  vector<double> nx = Linspace(20,120,6);

  // set problem information
  int Nt = 1000;
  double t0 = 0.0;
  double Tf = 1.0;
  double dtout = (Tf-t0)/Nt;
  double h = (Tf-t0)/Nt;

  // create ODE RHS function objects
  MyFlux flux;
  MySource source;
  
  // storage for error results
  vector<double> errs(nx.size());
  
 ///////// ERK3 /////////
  cout << "\nERK3 Method:\n";

  // loop over time step sizes
  for (int ih=0; ih<nx.size(); ih++) {

    // set up domain
    int a = -1;
    int b =  1;
    double dx = (b-a)/nx[ih];
    vector<double> x = Linspace(a + dx,b,nx[ih]);

    vector<double> y0;

    for(int i = 0; i < x.size(); i++)
      y0.push_back(sin(M_PI*x[i])); 

    MyRHS rhs(y0, flux, source, dx);
    ERK3Stepper ERK3(rhs, y0);
    

    // set the initial condition, initial time
    vector<double> y(y0);
    double tcur = t0;

    // reset maxerr
    double maxerr = 0.0;

    double stime = get_time();
    // loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt; istep++) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};
 
      // call the solver, update current time
      vector<double> tvals = ERK3.Evolve(tspan, h, y);
      tcur = tvals.back();     // last entry in tvals

      // calculate exact value at each time
      vector<double> exact = ytrue(tcur,x);
      // calculate the error
      vector<double> yerr = y - exact;
      double err = InfNorm(yerr);
      maxerr = std::max(maxerr, err);

    }
    double ftime = get_time();
    double rtime = ftime-stime;
    printf("testing time: %g\n", rtime);
    
    // compute convergence rate
    cout << "   n = " << nx[ih] << "\t  error = " << maxerr;
    errs[ih] = maxerr;
    if (ih > 0)
      cout << "\t  conv rate = " << -(log(errs[ih])-log(errs[ih-1]))/(log(nx[ih])-log(nx[ih-1]));
    cout << endl;

  }
    // output results to disk
    Matrix error(1, nx.size());
    Matrix n(1,nx.size());
    for(int j = 0; j < nx.size(); j++)
      n[j] = {nx[j]};
    for(int j = 0; j < nx.size(); j++)
      error[j] = {errs[j]};
    n.Write("nxvals.txt");
    error.Write("error.txt");


  return 0;
}
