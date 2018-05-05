#include "weno_mpi.hpp"


int ParallelDecomp::processLoc() {

  int ierr;
  if (nbE < 0) {
    int nbcoords[2];
    nbcoords[0] = pcoords[0]-1;
    nbcoords[1] = pcoords[1];
    ierr = MPI_Cart_rank(comm, nbcoords, &nbW);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
    nbcoords[0] = pcoords[0]+1;
    nbcoords[1] = pcoords[1];
    ierr = MPI_Cart_rank(comm, nbcoords, &nbE);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
    nbcoords[0] = pcoords[0];
    nbcoords[1] = pcoords[1]-1;
    ierr = MPI_Cart_rank(comm, nbcoords, &nbN);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
    nbcoords[0] = pcoords[0];
    nbcoords[1] = pcoords[1]+1;
    ierr = MPI_Cart_rank(comm, nbcoords, &nbS);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }  // end if (nbE < 0)
	return 0;
}

int ParallelDecomp::CommunicationWE(std::vector<std::vector<double>> &u) {
  
  // allocate send/recv arrays if NULL
  if (vrecvEp  == NULL)  vrecvEp  = new double[nyloc];
  if (vrecvEpp == NULL)  vrecvEpp = new double[nyloc];
  if (vsendEp  == NULL)  vsendEp  = new double[nyloc];
  if (vsendEpp == NULL)  vsendEpp = new double[nyloc];

  if (vrecvWp  == NULL)  vrecvWp  = new double[nyloc];
  if (vrecvWpp == NULL)  vrecvWpp = new double[nyloc];
  if (vsendWp  == NULL)  vsendWp  = new double[nyloc];
  if (vsendWpp == NULL)  vsendWpp = new double[nyloc]; 



  // initialize send/recv buffers
  for (int i=0; i<nyloc; i++)  vsendWp[i]  = u[i][0];
  for (int i=0; i<nyloc; i++)  vsendWpp[i] = u[i][1];
  for (int i=0; i<nyloc; i++)  vsendEp[i]  = u[i][nxloc-1];
  for (int i=0; i<nyloc; i++)  vsendEpp[i] = u[i][nxloc-2];
  for (int i=0; i<nyloc; i++)  vrecvEp[i]  = 0.0;
  for (int i=0; i<nyloc; i++)  vrecvEpp[i] = 0.0;
  for (int i=0; i<nyloc; i++)  vrecvWp[i]  = 0.0;
  for (int i=0; i<nyloc; i++)  vrecvWpp[i] = 0.0;

  int ierr;
  // open send/receive channels
  MPI_Request req[8];
  ierr = MPI_Irecv(vrecvWp, nyloc, MPI_DOUBLE, nbW, 2, comm, &(req[0]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvWpp, nyloc, MPI_DOUBLE, nbW, 4, comm, &(req[1]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvEp, nyloc, MPI_DOUBLE, nbE, 6, comm, &(req[2]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvEpp, nyloc, MPI_DOUBLE, nbE, 8, comm, &(req[3]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  ierr = MPI_Isend(vsendEp, nyloc, MPI_DOUBLE, nbE, 2, comm, &(req[4]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendEpp, nyloc, MPI_DOUBLE, nbE, 4, comm, &(req[5]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendWp, nyloc, MPI_DOUBLE, nbW, 6, comm, &(req[6]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendWpp, nyloc, MPI_DOUBLE, nbW, 8, comm, &(req[7]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // wait until all communications finish
  MPI_Status stat[8];
  ierr = MPI_Waitall(8, req, stat);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Waitall = %i\n",ierr);
    MPI_Abort(comm, 1);
  }


  return 0;

}

  
int ParallelDecomp::CommunicationNS(std::vector<std::vector<double>> &u) {

  
  if (vrecvNp  == NULL)  vrecvNp  = new double[nxloc];
  if (vrecvNpp == NULL)  vrecvNpp = new double[nxloc];
  if (vsendNp  == NULL)  vsendNp  = new double[nxloc];
  if (vsendNpp == NULL)  vsendNpp = new double[nxloc];

  if (vrecvSp  == NULL)  vrecvSp  = new double[nxloc];
  if (vrecvSpp == NULL)  vrecvSpp = new double[nxloc];
  if (vsendSp  == NULL)  vsendSp  = new double[nxloc];
  if (vsendSpp == NULL)  vsendSpp = new double[nxloc];


  // initialize send/recv buffers
  for (int i=0; i<nxloc; i++)  vsendSp[i]  = u[nyloc-1][i];
  for (int i=0; i<nxloc; i++)  vsendSpp[i] = u[nyloc-2][i];
  for (int i=0; i<nxloc; i++)  vsendNp[i]  = u[0][i];
  for (int i=0; i<nxloc; i++)  vsendNpp[i] = u[1][i];
  for (int i=0; i<nxloc; i++)  vrecvNp[i]  = 0.0;
  for (int i=0; i<nxloc; i++)  vrecvNpp[i] = 0.0;
  for (int i=0; i<nxloc; i++)  vrecvSp[i]  = 0.0;
  for (int i=0; i<nxloc; i++)  vrecvSpp[i] = 0.0;

  int ierr;
  // open send/receive channels
  MPI_Request req[8];
  ierr = MPI_Irecv(vrecvSp, nxloc, MPI_DOUBLE, nbS, 2, comm, &(req[0]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvSpp, nxloc, MPI_DOUBLE, nbS, 4, comm, &(req[1]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvNp, nxloc, MPI_DOUBLE, nbN, 6, comm, &(req[2]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvNpp, nxloc, MPI_DOUBLE, nbN, 8, comm, &(req[3]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }


  ierr = MPI_Isend(vsendNp, nxloc, MPI_DOUBLE, nbN, 2, comm, &(req[4]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendNpp, nxloc, MPI_DOUBLE, nbN, 4, comm, &(req[5]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendSp, nxloc, MPI_DOUBLE, nbS, 6, comm, &(req[6]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendSpp, nxloc, MPI_DOUBLE, nbS, 8, comm, &(req[7]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // wait until all communications finish
  MPI_Status stat[8];
  ierr = MPI_Waitall(8, req, stat);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Waitall = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  return 0;
}

int ParallelDecomp::CommunicationFluxWE(std::vector<std::vector<double>> & rf,
                      	std::vector<std::vector<double>> & lf){

  if (vsendEp == NULL)  vsendEp  = new double[nyloc];
  if (vsendWp == NULL)  vsendWp  = new double[nyloc];
  if (vrecvEp == NULL)  vrecvEp  = new double[nyloc];
  if (vrecvWp == NULL)  vrecvWp  = new double[nyloc];


  // initialize send/recv buffers
  for (int i=0; i<nyloc; i++)  vsendEp[i]  = rf[i][nxloc-1];
  for (int i=0; i<nyloc; i++)  vsendWp[i]  = lf[i][0];
  for (int i=0; i<nyloc; i++)  vrecvEp[i]  = 0.0;
  for (int i=0; i<nyloc; i++)  vrecvWp[i]  = 0.0;

  int ierr;
  // open send/receive channels
  MPI_Request req[4];
  ierr = MPI_Irecv(vrecvWp, nyloc, MPI_DOUBLE, nbW, 2, comm, &(req[0]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvEp, nyloc, MPI_DOUBLE, nbE, 4, comm, &(req[1]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  
  ierr = MPI_Isend(vsendEp, nyloc, MPI_DOUBLE, nbE, 2, comm, &(req[2]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendWp, nyloc, MPI_DOUBLE, nbW, 4, comm, &(req[3]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // wait until all communications finish
  MPI_Status stat[4];
  ierr = MPI_Waitall(4, req, stat);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Waitall = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  return 0;
 }


int ParallelDecomp::CommunicationFluxNS(std::vector<std::vector<double>> &rf,
                      	std::vector<std::vector<double>> & lf){

  if (vsendNp  == NULL)  vsendNp  = new double[nxloc];
  if (vsendSp  == NULL)  vsendSp  = new double[nxloc];
  if (vrecvNp  == NULL)  vrecvNp  = new double[nxloc];
  if (vrecvSp  == NULL)  vrecvSp  = new double[nxloc];

  // initialize send/recv buffers
  for (int i=0; i<nxloc; i++)  vsendNp[i] = lf[0][i];
  for (int i=0; i<nxloc; i++)  vsendSp[i] = rf[nyloc-1][i];
  for (int i=0; i<nxloc; i++)  vrecvNp[i] = 0.0;
  for (int i=0; i<nxloc; i++)  vrecvSp[i] = 0.0;


  int ierr;
  // open send/receive channels
  MPI_Request req[4];
  ierr = MPI_Irecv(vrecvNp, nxloc, MPI_DOUBLE, nbN, 2, comm, &(req[0]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Irecv(vrecvSp, nxloc, MPI_DOUBLE, nbS, 4, comm, &(req[1]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendSp, nxloc, MPI_DOUBLE, nbS, 2, comm, &(req[2]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Isend(vsendNp, nxloc, MPI_DOUBLE, nbN, 4, comm, &(req[3]));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Isend = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // wait until all communications finish
  MPI_Status stat[4];
  ierr = MPI_Waitall(4, req, stat);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Waitall = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
	

  return 0;


}







