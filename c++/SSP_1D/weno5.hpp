/* Fifth order WENO method file.

   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

#ifndef WENO5_DEFINED__
#define WENO5_DEFINED__

// Inclusions
#include <vector>       
#include "matrix.hpp"
#include "rhs.hpp"

// WENO5
class WEN05 {

 private:
  RHSFlux *flux;                                                 // pointer to Flux function
  RHSSource *source;                                             // pointer to source function
  std::vector<double> lf, rf, w, v, uf, s, f, df, vp, vn;  // reused vectors
  double dx, a;
  double umm,um,u,up,upp,ep,epsilon;
  double d0n,d1n,d2n,B0n,B1n,B2n,alpha0n,alpha1n,alpha2n,alphasumn,w0n,w1n,w2n,p0n,p1n,p2n;
  int N;
  
 public:
  // constructor (sets RHS function pointer, allocates local data)
  WEN05(std::vector<double>& u_, RHSFlux& flux_, RHSSource& source_, double dx_);


  int periodicBoundary(int, int, std::vector<double>&);
  // Reconstruction of the stenicls
  int reconstruction(int, std::vector<double>&, std::vector<double>&);

  // Compute Right Flux
  std::vector<double> residual(std::vector<double>&);

  int fluxSplitting(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);


};

#endif
