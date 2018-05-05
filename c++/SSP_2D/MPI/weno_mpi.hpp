/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */



// Parallelism class
#ifndef PARALLEL_DECOMP_DEFINED__
#define PARALLEL_DECOMP_DEFINED__


#include "mpi.h"
#include <vector>
#include <iostream>   

class ParallelDecomp{

 public:

  MPI_Comm comm;
  int pdims[2];
  int periodic[2];
  int pcoords[2];
  int nxloc;
  int nyloc;
  int myid;
  int numprocs;
  int nbE;
  int nbW;
  int nbN;
  int nbS;

  // For movement in x-direction
  double *vrecvEp;
  double *vrecvEpp;
  double *vsendEp;
  double *vsendEpp;
  double *vrecvWp;
  double *vrecvWpp;
  double *vsendWp;
  double *vsendWpp;

  // For movement in y-direction
  double *vrecvNp;
  double *vrecvNpp;
  double *vsendNp;
  double *vsendNpp;
  double *vrecvSp;
  double *vrecvSpp;
  double *vsendSp;
  double *vsendSpp;

  ParallelDecomp() {
    nbE = nbW = nbN = nbS = -1;
    vrecvEp = vrecvEpp = vsendEp = vsendEpp = NULL;
    vrecvNp = vrecvNpp = vsendNp = vsendNpp = NULL;
    vrecvSp = vrecvSpp = vsendSp = vsendSpp = NULL;
    vrecvWp = vrecvWpp = vsendWp = vsendWpp = NULL;
  }
  ~ParallelDecomp() {

    // Moving along the x direction

    // East
    if(vrecvEp  != NULL) delete[] vrecvEp;
    if(vrecvEpp != NULL) delete[] vrecvEpp;
    if(vsendEp  != NULL) delete[] vsendEp;
    if(vsendEpp != NULL) delete[] vsendEpp;

    // West
    if(vsendWp  != NULL) delete[] vsendWp;
    if(vsendWpp != NULL) delete[] vsendWpp;
    if(vrecvWp  != NULL) delete[] vrecvWp;
    if(vrecvWpp != NULL) delete[] vrecvWpp;

    // Moving along the y direction

    // North
    if(vrecvNp  != NULL) delete[] vrecvNp;
    if(vrecvNpp != NULL) delete[] vrecvNpp;
    if(vsendSp  != NULL) delete[] vsendSp;
    if(vsendSpp != NULL) delete[] vsendSpp;

    // South
    if(vrecvSp  != NULL) delete[] vrecvSp;
    if(vrecvSpp != NULL) delete[] vrecvSpp;
    if(vsendNp  != NULL) delete[] vsendNp;
    if(vsendNpp != NULL) delete[] vsendNpp;

  }

  int CommunicationWE(std::vector<std::vector<double>> &);
  int CommunicationNS(std::vector<std::vector<double>> &);
  int CommunicationFluxWE(std::vector<std::vector<double>> &,
                          std::vector<std::vector<double>> &);
  int CommunicationFluxNS(std::vector<std::vector<double>> &,
                          std::vector<std::vector<double>> &);
  int processLoc();
  
};  // end parallel_decomp


#endif