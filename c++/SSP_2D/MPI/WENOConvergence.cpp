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
#include "mpi.h"
#include "weno_mpi.hpp"

using namespace std;


// ODE RHS function
class MyRHS: public RHSFunction {
private:
  WEN05 *weno;
public:
  
  MyRHS(vector<vector<double>> &u_, RHSFlux& flux_, RHSSource& source_, double dx_,
       double dy_,ParallelDecomp& p){
    weno = new WEN05(u_, flux_, source_, dx_, dy_, p);
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
      for(int j = 0; j < y[i].size(); j++)
        f[i][j] = y[i][j];

    for(int i = 0; i < y.size(); i++)
      for(int j = 0; j < y[i].size(); j++)
        df[i][j] = 1;

    return 0;
  }
};

// ODE Source function
class MySource : public RHSSource{
public:
  int Evaluate(vector<vector<double>>& y, vector<vector<double>>& f) {
    for(int i = 0; i  < y.size(); i++)
      for(int j = 0; j < y[i].size(); j++)
        f[i][j] = 0;
    return 0;
  }
};

// Convenience function for analytical solution
vector<vector<double>> ytrue(const double t, vector<double>& x, vector<double>& y) { 
  vector<vector<double>> yt;
  for(int i = 0; i  < y.size(); i++){
    yt.push_back(vector<double>());
    for(int j = 0; j  < x.size(); j++)
      yt[i].push_back(sin(2.0*M_PI*(x[j] - t))*sin(2.0*M_PI*(y[i] - t)));
  }
  return yt;
};

// Function prototypes
void testConvergence(vector<int>); 
void showScalability(vector<int>, vector<int>);
void getInfo();

// main routine
int main(int argc, char *argv[] ) {

  // Define local variables
  int ierr = 0, numprocs,myid;
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Init = " << ierr << "\n";
    return 1;
  }

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_size = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_rank = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  
  if(argc > 1){
    vector<int> contents; 
    string input;
    for (int i = 1; i < argc; i++){
      
      if (isdigit(argv[i][0])){
        contents.push_back(atoi(argv[i]));
      }
      else
        input = argv[i];
    }

    if(input == "c" || input == "C"){
      vector<int> P = { contents[0], contents[1] };
      testConvergence(P);
    }
    else if(input == "s" || input == "S"){
      vector<int> N = { contents[0], contents[1] };
      vector<int> P = { contents[2], contents[3] };
      showScalability(N,P);
    }
    else if(input == "i" || input == "I"){
      getInfo();
    }
  }
  
    // finalize MPI
  ierr = MPI_Finalize();

  return 0;
}
void getInfo(){
 #ifdef _OPENMP
  printf("Number of Processors = %i\n",omp_get_num_procs() );
  printf("Number of Threads = %i\n", omp_get_num_threads());
  printf("Max Threads = %i\n\n", omp_get_max_threads());
 #else
  int ierr = 0, numprocs, myid;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_size = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_rank = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(myid == 0)
    printf("Number of Processors = %i\n", numprocs);
 #endif
}
void showScalability(vector<int> N, vector<int> P){

  int i, j, ih, istep, Nt, a, b, px, py, ierr;
  double t0, Tf, dtout, h, dx, dy, tcur, maxerr, allmaxerr;
  double stime, ftime, rtime, err, nx, ny, cfl, maxflux, allmaxflux;
  vector<double> x, y, xl, yl, tspan, tvals;
  vector<vector<double>> u, u0, allu0, exact, yerr;
  ParallelDecomp p2d;

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(p2d.numprocs));
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_size = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(p2d.myid));
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_rank = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Number of Mesh Pts 
  nx = N[0];
  ny = N[1];
  // Number of processes in each direction
  px = P[0];
  py = P[1];

  // create ODE RHS function objects
  MyFlux flux;
  MySource source;

  // check for correct process dimension
  if (p2d.numprocs != px*py) {
    std::cerr << " Incorrect processor layout, px = " << px << ", py = " << py 
            << ", numprocs = " << p2d.numprocs << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  // perform parallel decomposition
  p2d.pdims[0] = px;
  p2d.pdims[1] = py;
  p2d.periodic[0] = 1;
  p2d.periodic[1] = 1;
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, p2d.pdims, p2d.periodic, 0, &(p2d.comm));
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Cart_create = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Cart_get(p2d.comm, 2, p2d.pdims, p2d.periodic, p2d.pcoords);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Cart_get = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int is, ie, js, je;
  is = p2d.pcoords[0]*nx/px;
  ie = (p2d.pcoords[0]+1)*nx/px-1;
  js = p2d.pcoords[1]*ny/py;
  je = (p2d.pcoords[1]+1)*ny/py-1;
  
  // allocate arrays
  stime = MPI_Wtime();
  p2d.nxloc = ie-is+1;
  p2d.nyloc = je-js+1;

  // Find process locations
  p2d.processLoc();

    ///////// ERK3 /////////
  if(p2d.myid == 0)
    cout << "ERK3 Method:\n";

  a =  0;
  b =  1;
  dx = (b-a)/nx;
  dy = (b-a)/ny;

  // // set mesh points
  for (i=is; i<=ie; i++) {
    x.push_back(a + dx*(i+1));
  }
  for (j=js; j<=je; j++) {
    y.push_back(a + dy*(j+1));
  }

  u0 = {};
  for(i = 0; i < p2d.nyloc; i++){
    u0.push_back(vector<double>());
    for(j = 0; j < p2d.nxloc; j++)
      u0[i].push_back(sin(2.0*M_PI*x[j])*sin(2.0*M_PI*y[i])); 
  }


  vector<vector<double>> w = u0;
  vector<vector<double>> f = u0;
  vector<vector<double>> df = u0;

  // set problem information
  if (flux.Evaluate(w,f,df) != 0) {
    std::cerr << "Step: Error in Flux function\n";
  }

  for(i = 0; i < p2d.nyloc; i++)
    for(j = 0; j < p2d.nxloc; j++)
        maxflux = std::max(maxflux,df[i][j]);

  cfl = 0.4;
  t0 = 0.0;
  Tf = 1.0;
  Nt = (int) (Tf-t0)*maxflux/(cfl*dx);
  dtout = (Tf-t0)/Nt;
  h = (Tf-t0)/Nt;

  MyRHS * rhs = new MyRHS(u0, flux, source, dx, dy, p2d);
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
    for(i = 0; i < p2d.nyloc; i++){
      yerr.push_back(vector<double>());
      for(j = 0; j < p2d.nxloc; j++)
        yerr[i].push_back(u[i][j] - exact[i][j]);
    }

    // infinity norm
    for(i = 0; i < p2d.nyloc; i++)
      for(j = 0; j < p2d.nxloc; j++)
        err = std::max(err,yerr[i][j]);
 
    maxerr = std::max(maxerr, err);

  }

  ierr = MPI_Allreduce(&maxerr, &allmaxerr, 1, MPI_DOUBLE, MPI_MAX, p2d.comm);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Allreduce = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // End time
  ftime = get_time();
  rtime = ftime-stime;

  // compute convergence rate
  if(p2d.myid==0){
    printf("\tnx = %g\n", nx);
    printf("\tny = %g\n", ny);
    printf("\tNt = %i\n", Nt);
    printf("\tlocal nx = %i\n", p2d.nxloc);
    printf("\tlocal ny = %i\n", p2d.nyloc);
    printf("\tlocal error = %g\n", maxerr);
    printf("\tactual error = %g\n", allmaxerr);
    printf("\ttesting time: %g\n", rtime);
  }
  

  delete rhs;
  delete ERK3;

}

void testConvergence(vector<int> P){

  int i, j, ih, istep, Nt, a, b, px, py, ierr;
  double t0, Tf, dtout, h, dx, dy, tcur, maxerr, allmaxerr;
  double stime, ftime, rtime, err, cfl, maxflux, allmaxflux;
  vector<double> x, y, xl, yl, tspan, tvals,nx;
  vector<vector<double>> u, u0, allu0, exact, yerr;
  ParallelDecomp p2d;

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(p2d.numprocs));
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_size = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(p2d.myid));
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_rank = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Number of processes in each direction
  px = P[0];
  py = P[1];


  if (p2d.numprocs != px*py) {
    std::cerr << " Incorrect processor layout, px = " << px << ", py = " << py 
            << ", numprocs = " << p2d.numprocs << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  // perform parallel decomposition
  p2d.pdims[0] = px;
  p2d.pdims[1] = py;
  p2d.periodic[0] = 1;
  p2d.periodic[1] = 1;
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, p2d.pdims, p2d.periodic, 0, &(p2d.comm));
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Cart_create = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Cart_get(p2d.comm, 2, p2d.pdims, p2d.periodic, p2d.pcoords);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Cart_get = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Find process locations
  p2d.processLoc();

  // storage for error results
  vector<double> errs = {0, 0, 0, 0, 0};
  vector<double> nxs = {20, 40, 80, 160, 320};
  vector<double> nys = {20, 40, 80, 160, 320};

  ///////// ERK3 /////////
  if(p2d.myid == 0)
    cout << "ERK3 Method:\n";

  for(int ih=0;ih<nxs.size();ih++){

    // Making processor grid
    int is, ie, js, je;
    is = p2d.pcoords[0]*nxs[ih]/px;
    ie = (p2d.pcoords[0]+1)*nxs[ih]/px-1;
    js = p2d.pcoords[1]*nys[ih]/py;
    je = (p2d.pcoords[1]+1)*nys[ih]/py-1;

    // allocate arrays
    stime = MPI_Wtime();
    p2d.nxloc = ie-is+1;
    p2d.nyloc = je-js+1;


    // set problem information
    Nt = 5000;
    t0 = 0.0;
    Tf = 1.0;

    dtout = (Tf-t0)/Nt;
    h = (Tf-t0)/Nt;
    a =  0;
    b =  1;
    dx = (b-a)/nxs[ih];
    dy = (b-a)/nys[ih];

    // // set mesh points
    x = {};
    for (i=is; i<=ie; i++) {
      x.push_back(a + dx*(i+1));
    }
    y = {};
    for (j=js; j<=je; j++) {
      y.push_back(a + dy*(j+1));
    }

    u0 = {};
    for(i = 0; i < p2d.nyloc; i++){
      u0.push_back(vector<double>());
      for(j = 0; j < p2d.nxloc; j++)
        u0[i].push_back(sin(2.0*M_PI*x[j])*sin(2.0*M_PI*y[i])); 
    }



        // create ODE RHS function objects
    MyFlux flux;
    MySource source;
    MyRHS * rhs = new MyRHS(u0, flux, source, dx, dy, p2d);
    ERK3Stepper *ERK3 = new ERK3Stepper(rhs, u0);  

    // set the initial condition, initial time
    u = u0;
    tcur = t0;

    // reset maxerr
    err = 0;
    maxerr = 0.0;  
    tvals = {};
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
        for(i = 0; i < p2d.nyloc; i++){
          yerr.push_back(vector<double>());
          for(j = 0; j < p2d.nxloc; j++)
            yerr[i].push_back(u[i][j] - exact[i][j]);
        }

        // infinity norm
        for(i = 0; i < p2d.nyloc; i++)
          for(j = 0; j < p2d.nxloc; j++)
            err = std::max(err,yerr[i][j]);
     
        maxerr = std::max(maxerr, err);

      }

    ierr = MPI_Allreduce(&maxerr, &allmaxerr, 1, MPI_DOUBLE, MPI_MAX, p2d.comm);
    if (ierr != MPI_SUCCESS) {
      std::cerr << " error in MPI_Allreduce = " << ierr << "\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    maxerr = allmaxerr;

    ftime = get_time();
    rtime = ftime-stime;
    if(p2d.myid == 0){
      
      // compute convergence rate
      printf("\tnx = %g", nxs[ih]);
      printf(",\tny = %g", nys[ih]);
      printf(",\terror = %g", maxerr);

      errs[ih] = maxerr;
      if (ih > 0)
        printf("\t  conv rate = %g", -(log(errs[ih])-log(errs[ih-1]))/(log(nxs[ih])-log(nxs[ih-1])));
      printf("\n");
    }

    delete rhs;
    delete ERK3;

  } 


}

