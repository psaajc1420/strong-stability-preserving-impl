/* Main routine to test the stability of strong stability Runge-Kutta 
   time discretizations of order 2, 3, and 4 against the corresponding 
   order of explicit Runge-Kutta methods 

      
     y' = L(u), t in [0,1],         // ODE
     L(u) = u^2/(u^2 + (1 - u)^2)   // Buckey-Leverett Flux
     y(0) = {                       // Initial Condition
              2 ,  x <= x_mid
              1 ,  else 
            }

   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <iomanip>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "erk2.hpp"
#include "erk3.hpp"
#include "erk4.hpp"
#include "ssp2.hpp"
#include "ssp3.hpp"
#include "ssp4.hpp"
#include "weno5.hpp"


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
    f = weno->reconstruction(y);
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
      f[i] = pow(y[i],2.0)/(pow(y[i],2.0) + pow((1.0 - y[i]),2));

    for(int i = 0; i < y.size(); i++)
      df[i] = -(2.0*(y[i] -1.0)*y[i])/pow((2.0*pow(y[i],2) - 2.0*y[i] + 1.0),2);

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


// main routine
int main() {

  // create ODE RHS function objects
  MyFlux flux;
  MySource source;

   // set up domain
  double a = -1;
  double b =  1;
  double nx = 160;
  double dx = (b-a)/nx;
  vector<double> x = Linspace(a + dx,b,nx);

  double t0 = 0.0;
  double Tf = 1.0;

  vector<double> y0;

  double mid = 0.5*(a + b);
  for(int i = 0; i < x.size(); i++){
    if(x[i] <= mid)
    y0.push_back(2.0); 
    else
    y0.push_back(1.0);
  }

  MyRHS rhs(y0, flux, source, dx);
  ERK2Stepper ERK2(rhs, y0);
  ERK3Stepper ERK3(rhs, y0);
  ERK4Stepper ERK4(rhs, y0);
  SSP2Stepper SSP2(rhs, y0);
  SSP3Stepper SSP3(rhs, y0);
  SSP4Stepper SSP4(rhs, y0);

  // set the initial condition, initial time
  vector<double> y;
  vector<double> f(y0);
  vector<double> df(y0);

  vector<double> CFLS2 = {1.0, 1.5, 0.5};
  vector<double> CFLS3 = {1.0, 1.5, 0.5};
  vector<double> CFLS4 = {2.0, 1.4, 1.0};

  // set problem information
  if (flux.Evaluate(y,f,df) != 0) {
    std::cerr << "Step: Error in Flux function\n";
  }

  double maxflux = max(abs(df));

    ///////// ERK2 /////////
      cout << "\nERK2 Methods:\n";
  for(int ii = 0; ii < CFLS2.size(); ii++){
      
      double CFL = CFLS2[ii];
      double tcur = t0;
      double dtout = CFL*dx/maxflux;
      double h = dtout;
      double Nt = (Tf-t0)/h;
      Matrix ycal(Nt+1,nx);
      y = y0;
      for(int j = 0; j < y.size(); j++ )
        ycal(0,j) = y[j];
      

      // loop over output step sizes: call solver and output error
      for (int istep=0; istep<Nt; istep++) {

        // set the time interval for this solve
        vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

        // call the solver, update current time
        vector<double> tvals = ERK2.Evolve(tspan, h, y);
        tcur = tvals.back();     // last entry in tvals

        // store approximated solution and exact solution for each time
        for(int jstep = 0; jstep < y.size(); jstep++)
          ycal(istep +1,jstep) = y[jstep];    

      }

      //output results to disk
      if(ii == 0){
        Matrix xt(1, x.size());
        Matrix yo(1, y0.size());
        for(int j = 0; j < x.size(); j++)
          xt[j] = {x[j]};
        for(int j = 0; j < y0.size(); j++)
          yo[j] = {y0[j]};

        xt.Write("x_Nonlinear.txt");
        yo.Write("y0_Nonlinear.txt");
      }
      
      char fname[35];
      sprintf(fname,"ycalc_NonlinearBL_erk2_CFL(%g).txt",CFL);
      ycal.Write(fname);
    }


         ///////// SSP2 /////////
      cout << "\nSSP2 Method:\n";

      // reset the initial condition, initial time
    for(int ii = 0; ii < CFLS2.size(); ii++){
      double CFL = CFLS2[ii];
      double tcur = t0;
      double dtout = CFL*dx/maxflux;
      double h = dtout;
      double Nt = (Tf-t0)/h;     
      Matrix ycal(Nt+1,nx);
      y = y0;

      for(int j = 0; j < y.size(); j++ )
        ycal(0,j) = y[j];


      // loop over output step sizes: call solver and output error
      for (int istep=0; istep<Nt; istep++) {

        // set the time interval for this solve
        vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

        // call the solver, update current time
        vector<double> tvals = SSP2.Evolve(tspan, h, y);
        tcur = tvals.back();     // last entry in tvals

        // store approximated solution and exact solution for each time
        for(int jstep = 0; jstep < y.size(); jstep++)
          ycal(istep +1,jstep) = y[jstep];
        

      }

      char fname[35];
      sprintf(fname,"ycalc_NonlinearBL_ssp2_CFL(%g).txt",CFL);
      ycal.Write(fname);
    }


      ///////// ERK3 /////////
      cout << "\nERK3 Method:\n";

    for(int ii = 0; ii < CFLS3.size(); ii++){

        double CFL = CFLS3[ii];
        double tcur = t0;
        double dtout = CFL*dx/maxflux;
        double h = dtout;
        double Nt = (Tf-t0)/h;
        Matrix ycal(Nt+1,nx);
        y = y0;

        for(int j = 0; j < y.size(); j++ )
          ycal(0,j) = y[j];
      

      // loop over output step sizes: call solver and output error
      for (int istep=0; istep<Nt; istep++) {

        // set the time interval for this solve
        vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

        // call the solver, update current time
        vector<double> tvals = ERK3.Evolve(tspan, h, y);
        tcur = tvals.back();     // last entry in tvals

        // store approximated solution and exact solution for each time
        for(int jstep = 0; jstep < y.size(); jstep++)
          ycal(istep +1,jstep) = y[jstep];
        

      }
      char fname[35];
      sprintf(fname,"ycalc_NonlinearBL_erk3_CFL(%g).txt",CFL);
      ycal.Write(fname);
    }


    ///////// SSP3 /////////
    cout << "\nSSP3 Method:\n";

    for(int ii = 0; ii < CFLS3.size(); ii++){

      double CFL = CFLS3[ii];
      double tcur = t0;
      double dtout = CFL*dx/maxflux;
      double h = dtout;
      double Nt = (Tf-t0)/h;
      Matrix ycal(Nt+1,nx);
      y = y0;
      for(int j = 0; j < y.size(); j++ )
        ycal(0,j) = y[j];


    // loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt; istep++) {

      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver, update current time
      vector<double> tvals = SSP3.Evolve(tspan, h, y);
      tcur = tvals.back();     // last entry in tvals

      // store approximated solution and exact solution for each time
      for(int jstep = 0; jstep < y.size(); jstep++)
        ycal(istep +1,jstep) = y[jstep];
    

    }

    char fname[35];
    sprintf(fname,"ycalc_NonlinearBL_ssp3_CFL(%g).txt",CFL);
    ycal.Write(fname);


  }


  ///////// ERK4 /////////
  cout << "\nERK4 Method:\n";

  for(int ii = 0; ii < CFLS4.size(); ii++){

      double CFL = CFLS4[ii];
      double tcur = t0;
      double dtout = CFL*dx/maxflux;
      double h = dtout;
      double Nt = (Tf-t0)/h;
      y = y0;

      Matrix ycal(Nt+1,nx);

      for(int j = 0; j < y.size(); j++ )
        ycal(0,j) = y[j];


      // loop over output step sizes: call solver and output error
      for (int istep=0; istep<Nt; istep++) {

        // set the time interval for this solve
        vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

        // call the solver, update current time
        vector<double> tvals = ERK4.Evolve(tspan, h, y);
        tcur = tvals.back();     // last entry in tvals

        // store approximated solution and exact solution for each time
        for(int jstep = 0; jstep < y.size(); jstep++)
          ycal(istep +1,jstep) = y[jstep];
      

      }
       char fname[35];
       sprintf(fname,"ycalc_NonlinearBL_erk4_CFL(%g).txt",CFL);
       ycal.Write(fname);
    }

    ///////// SSP4 /////////
    cout << "\nSSP4 Method:\n";
    
    for(int ii = 0; ii < CFLS4.size(); ii++){

      double CFL = CFLS4[ii];
      double tcur = t0;
      double dtout = CFL*dx/maxflux;
      double h = dtout;
      double Nt = (Tf-t0)/h;
      Matrix ycal(Nt+1,nx);
      y = y0;
    for(int j = 0; j < y.size(); j++ )
      ycal(0,j) = y[j];
    

    // loop over output step sizes: call solver and output error
    for (int istep=0; istep<Nt; istep++) {

      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver, update current time
      vector<double> tvals = SSP4.Evolve(tspan, h, y);
      tcur = tvals.back();     // last entry in tvals

      // store approximated solution and exact solution for each time
      for(int jstep = 0; jstep < y.size(); jstep++)
        ycal(istep +1,jstep) = y[jstep];
      

    }
    char fname[35];
    sprintf(fname,"ycalc_NonlinearBL_ssp4_CFL(%g).txt",CFL);
    ycal.Write(fname);
    
  }    
  

  return 0;
}
