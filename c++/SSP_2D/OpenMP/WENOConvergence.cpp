/* Main routine to test convergence of fifth order WENO method with 
   the a linear flux

     y' = L(u), t in [0,1],       // ODE
     L(u) = a                     // Linear Flux
     y(0) = sin(pi*x)*sin(pi*y)   // Initial Condition
   
   Jacob Cadena
   Math 6370 @ SMU
   Fall 2016  */

#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "matrix.hpp"
#include "rhs.hpp"
#include "erk3.hpp"
#include "weno5.hpp"
#include "get_time.h"
#include <omp.h> 


using namespace std;


// ODE RHS function
class MyRHS: public RHSFunction {
private:
  WEN05 *weno;
public:
  
  MyRHS(vector<vector<double>> &u_, RHSFlux& flux_, RHSSource& source_, double dx_, double dy_){
    weno = new WEN05(u_, flux_, source_, dx_, dy_);
  }
  int Evaluate(double t, vector<vector<double>> &y, vector<vector<double>>& f) { 
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
  int Evaluate(vector<vector<double>>& y,vector<vector<double>>& f, vector<vector<double>>& df) {
    for(int i = 0; i < y.size(); i++)
      for(int j = 0; j  < y[i].size(); j++)
        f[i][j] = y[i][j];

    for(int i = 0; i < y.size(); i++)
      for(int j = 0; j  < y[i].size(); j++)
        df[i][j] = 1;

    return 0;
  }
};

// ODE Source function
class MySource : public RHSSource{
public:
  int Evaluate(vector<vector<double>>& y, vector<vector<double>>& f) {
    for(int i = 0; i  < y.size(); i++)
      for(int j = 0; j  < y[i].size(); j++)
        f[i][j] = 0;
    return 0;
  }
};

// Convenience function for analytical solution
vector<vector<double>> ytrue(const double t, vector<double>& x, vector<double>& y) { 
  vector<vector<double>> yt;
  for(int i = 0; i  < x.size(); i++){
    yt.push_back(vector<double>());
    for(int j = 0; j  < y.size(); j++)
      yt[i].push_back(sin(2.0*M_PI*(x[i] - t))*sin(2.0*M_PI*(y[j] - t)));
  }
  return yt;
};


void testConvergence();
void showScalability(int,int);
void getInfo();

// main routine
int main(int argc, char *argv[] ) {

  if(argc > 1){
    vector<int> N; 
    for (int i = 1; i < argc; i++){
      string input;
      if (isdigit(argv[i][0])){
        N.push_back(atoi(argv[i]));
      }
      else
        input = argv[i];

      if(input == "c" || input == "C"){
        testConvergence();
      }
      else if(input == "s" || input == "S"){
        showScalability(N[0], N[1]);
      }
      else if(input == "i" || input == "I"){
        getInfo();
      }
    }
  }
  else{

    while(true){
      string input = "";
      printf("Enter (c) to show Convergence\n");
      printf("Enter (s) to show Scalability\n");
      printf("Enter (i) to get Info \n");
      printf("Enter (e) to Exit \n");

      getline(cin, input);

      printf("\n\n");
      if(input == "c" || input == "C"){
        testConvergence();
        break;
      }
      else if(input == "s" || input == "S"){
        printf("Enter nx:  ");
        getline(cin, input);
        int nx = stoi(input);
        printf("\nEnter ny:  ");
        getline(cin, input);
        int ny = stoi(input);
        showScalability(nx,ny);
        break;
      }
      else if(input == "i" || input == "I"){
        getInfo();
      }
      else if(input == "E" || input == "e"){
        printf("Goodbye :) \n");
        break;
      }
      else{
        printf("Invalid Input! Try Again :) \n");
      }
    }
  }
  


  return 0;
}
void getInfo(){
  printf("Number of Processors = %i\n",omp_get_num_procs() );
  printf("Number of Threads = %i\n", omp_get_num_threads());
  printf("Max Threads = %i\n\n", omp_get_max_threads());
}
void showScalability(int nx_, int ny_){

  int i, j, ih, istep, Nt, a, b;
  double t0, Tf, dtout, h, dx, dy, tcur, maxerr, stime, ftime, rtime, err, nx, ny, cfl, maxflux;
  vector<double> x, y, tspan, tvals;
  vector<vector<double>> u, u0, exact, yerr;

  // Number of Mesh Pts
  nx = nx_;
  ny = ny_;

  // create ODE RHS function objects
  MyFlux flux;
  MySource source;

  // storage for error results

  ///////// ERK3 /////////
  cout << "\nERK3 Method:\n";

  // set up domain
  a =  0;
  b =  1;
  dx = (b-a)/nx;
  dy = (b-a)/ny;
  x = Linspace(a + dx,b,nx);
  y = Linspace(a + dy,b,ny);

  u0 = {};
  for(i = 0; i < x.size(); i++){
    u0.push_back(vector<double>());
    for(j = 0; j < y.size(); j++)
      u0[i].push_back(sin(2.0*M_PI*x[i])*sin(2.0*M_PI*y[j])); 
  }

  for(i = 0; i < x.size(); i++){
    for(j = 0; j < y.size(); j++)
      cout << u0[i][j] << "  ";
    cout << "\n";
  }

  vector<vector<double>> w = u0;;
  vector<vector<double>> f = w;;
  vector<vector<double>> df = w;

  // set problem information
  if (flux.Evaluate(w,f,df) != 0) {
    std::cerr << "Step: Error in Flux function\n";
  }

  maxflux = InfNorm(df);

  cfl = 0.4;
  t0 = 0.0;
  Tf = 1.0;
  // Nt = 5000;
  Nt = (int) (Tf-t0)*maxflux/(cfl*dx);
  dtout = (Tf-t0)/Nt;
  h = (Tf-t0)/Nt;



  MyRHS * rhs = new MyRHS(u0, flux, source, dx, dy);
  ERK3Stepper *ERK3 = new ERK3Stepper(rhs, u0);  

  // set the initial condition, initial time
  u = u0;
  tcur = t0;

  // reset maxerr
  maxerr = 0.0;

  stime = get_time();    

  // loop over output step sizes: call solver and output error

    for (istep=0; istep<Nt; istep++){

      // set the time interval for this solve
      tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver, update current time
      tvals = ERK3->Evolve(tspan, h, u);

      tcur = tvals.back();     // last entry in tvals

      // calculate exact value at each time
      exact = ytrue(tcur, x, y);

      // calculate the error
      yerr = {};

      for(i = 0; i < x.size(); i++){
        yerr.push_back(vector<double>());
        for(j = 0; j < y.size(); j++)
          yerr[i].push_back(u[i][j] - exact[i][j]);
      }

      err = InfNorm(yerr);
      maxerr = std::max(maxerr, err);

    }

    // compute convergence rate
  printf("\tnx = %g\n", nx);
  printf("\tny = %g\n", ny);
  printf("\terror = %g\n", maxerr);

  ftime = get_time();
  rtime = ftime-stime;
  printf("\ttesting time: %g\n", rtime);
  

  delete rhs;
  delete ERK3;

}

void testConvergence(){

  int i, j, ih, istep, Nt, a, b;
  double t0, Tf, dtout, h, dx, dy, tcur, maxerr, stime, ftime, rtime, err;
  vector<double> nx, ny, errs, x, y, tspan, tvals;
  vector<vector<double>> u, u0, exact, yerr;

  nx = {20, 40, 60, 80, 100};
  ny = {20 ,40, 60, 80, 100};


  // set problem information
  Nt = 5000;
  t0 = 0.0;
  Tf = 1.0;

  dtout = (Tf-t0)/Nt;
  h = (Tf-t0)/Nt;

  // create ODE RHS function objects
  MyFlux flux;
  MySource source;

  // storage for error results
  errs = {0, 0, 0, 0, 0};

  ///////// ERK3 /////////
  cout << "\nERK3 Method:\n";
  
  // loop over grid sizes 
  for (ih=0; ih<nx.size(); ih++) {

    // set up domain
    a =  0;
    b =  1;
    dx = (b-a)/nx[ih];
    dy = (b-a)/ny[ih];
    x = Linspace(a + dx,b,nx[ih]);
    y = Linspace(a + dy,b,ny[ih]);

    u0 = {};
    for(i = 0; i < x.size(); i++){
      u0.push_back(vector<double>());
      for(j = 0; j < y.size(); j++)
        u0[i].push_back(sin(2.0*M_PI*x[i])*sin(2.0*M_PI*y[j])); 
    }

    MyRHS * rhs = new MyRHS(u0, flux, source, dx, dy);
    ERK3Stepper *ERK3 = new ERK3Stepper(rhs, u0);
    
    // set the initial condition, initial time
    u = u0;
    tcur = t0;

    // reset maxerr
    maxerr = 0.0;  

    stime = get_time();   

    // loop over output step sizes: call solver and output error

      for (istep=0; istep<Nt; istep++){

        // set the time interval for this solve
        tspan = {tcur, std::min(tcur + dtout, Tf)};

        // call the solver, update current time
        tvals = ERK3->Evolve(tspan, h, u);

        tcur = tvals.back();     // last entry in tvals

        // calculate exact value at each time
        exact = ytrue(tcur, x, y);

        // calculate the error
        yerr = {};

        for(i = 0; i < x.size(); i++){
          yerr.push_back(vector<double>());
          for(j = 0; j < y.size(); j++)
            yerr[i].push_back(u[i][j] - exact[i][j]);
        }

        err = InfNorm(yerr);
        maxerr = std::max(maxerr, err);

      }

    ftime = get_time();
    rtime = ftime-stime;
    printf("\ttesting time: %g\n", rtime);
    
    // compute convergence rate
    printf("\tnx = %g", nx[ih]);
    printf(",\tny = %g", ny[ih]);
    printf(",\terror = %g", maxerr);

    errs[ih] = maxerr;
    if (ih > 0)
      printf("\t  conv rate = %g", -(log(errs[ih])-log(errs[ih-1]))/(log(nx[ih])-log(nx[ih-1])));
    printf("\n");


    delete rhs;
    delete ERK3;

  }

}

