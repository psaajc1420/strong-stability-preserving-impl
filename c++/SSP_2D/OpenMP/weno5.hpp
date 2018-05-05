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
#include <omp.h> 

// WENO5
class WEN05 {

 private:

  RHSFlux *flux;                                                 // pointer to Flux function
  RHSSource *source;                                             // pointer to source function
  double dx,dy,a;
  int N, Nx, Ny;


  
 public:

  // Constructor (sets RHS function pointer, allocates local data) for 2D problem
  WEN05(std::vector<std::vector<double>> &u_, RHSFlux& flux_, RHSSource& source_, double dx_, double dy_);

  int reconstructionX(int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
  int reconstructionY(int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);


  // Compute 2D Residual
  std::vector<std::vector<double>> residual(std::vector<std::vector<double>>&);

  // Compute fluxsplitting
  void fluxSplitting(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&,
                     std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);


};

#endif