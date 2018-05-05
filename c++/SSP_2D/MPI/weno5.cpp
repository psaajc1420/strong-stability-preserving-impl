/* Fifth order WENO method file.

   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

// Inclusions    
#include "weno5.hpp"



WEN05::WEN05(std::vector<std::vector<double>> &u_, RHSFlux& flux_, RHSSource& source_
            ,double dx_, double dy_, ParallelDecomp& p) { 

  flux = &flux_;
  source = &source_;
  dx = dx_;
  dy = dy_;
  Ny = u_.size();
  Nx = u_[0].size();
  p2d = p;


}


void WEN05::fluxSplitting(std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &f, std::vector<std::vector<double>> &df,
                          std::vector<std::vector<double>> &up, std::vector<std::vector<double>> &un){

    int i, j;
    double a;

    // infinity norm
    for(i = 0; i < p2d.nyloc; i++)
      for(j = 0; j < p2d.nxloc; j++)
          a = std::max(a,df[i][j]);
  
    // # pragma omp parallel for collapse(2) private(i,j)
    for(i= 0; i < Ny; i++){
      for(j = 0; j < Nx; j++)
        up[i][j] = 0.5*(f[i][j] + a*u[i][j]);
    }

    // # pragma omp parallel for collapse(2) private(i,j)
    for(i = 0; i < Ny; i++){
      for(j = 0; j < Nx; j++)
        un[i][j] = 0.5*(f[i][j] - a*u[i][j]);
    }

}


std::vector<std::vector<double>> WEN05::residual(std::vector<std::vector<double>> &res) {

  int i, j;
  std::vector<std::vector<double>> w, vn, vp, f, df, dux, duy, rfx, rfy, lfx, lfy, s;

  w  = res;
  vn = res;
  vp = res;
  f  = res;
  df = res;
  dux = res;
  duy = res;
  rfx = res;
  rfy = res;
  lfx = res;
  lfy = res;
  s = res;

  // Calculate flux and dflux

  if (flux->Evaluate(w,f,df)!= 0) {
    std::cerr << "Step: Error in Flux function\n";
  }

  // Calculate Source 
  if(source->Evaluate(w,s)!=0) {
    std::cerr << "Step: Error in Source function\n";
  }
  // Along x
  // Lax-Friedrichs Flux Splitting
  fluxSplitting(w, f, df, vp, vn); 


  p2d.CommunicationWE(vp);
  reconstructionX( 1, vp, rfx);

  p2d.CommunicationWE(vn);
  reconstructionX(-1, vn, lfx);

  p2d.CommunicationFluxWE(rfx,lfx);
  // # pragma omp parallel for private(i,j)
  for(i = 0; i < Ny; i++){
    dux[i][0] = -(lfx[i][0] - lfx[i][1] + rfx[i][0] - p2d.vrecvWp[i])/dx + s[i][0];
    for(j = 1; j < Nx-1; j++)
      dux[i][j] = -(lfx[i][j] - lfx[i][j+1] + rfx[i][j] - rfx[i][j-1])/dx + s[i][j];
    dux[i][Nx-1] = -(lfx[i][Nx-1] - p2d.vrecvEp[i] + rfx[i][Nx-1] - rfx[i][Nx-2])/dx + s[i][Nx-1];
  }


  // Along y
  // Lax-Friedrichs Flux Splitting
  fluxSplitting(w, f, df, vp, vn); 

  p2d.CommunicationNS(vp);
  reconstructionY( 1, vp, rfy);

  p2d.CommunicationNS(vn);
  reconstructionY(-1, vn, lfy);


  p2d.CommunicationFluxNS(rfy,lfy);
  // # pragma omp parallel for private(i,j)
  for(j = 0; j < Nx; j++){
    duy[0][j] =  -(lfy[0][j] - lfy[1][j] + rfy[0][j] - p2d.vrecvNp[j])/dy + s[0][j];
    for(i =1; i < Ny-1; i++)
      duy[i][j] =  -(lfy[i][j] - lfy[i+1][j] + rfy[i][j] - rfy[i-1][j])/dy + s[i][j];
    duy[Ny-1][j] =  -(lfy[Ny-1][j] - p2d.vrecvSp[j] + rfy[Ny-1][j] - rfy[Ny-2][j])/dy + s[Ny-1][j];
  }


  // # pragma omp parallel for collapse(2) private(i,j)
  for(i = 0; i < Ny; i++)
    for(j = 0; j < Nx; j++)
      res[i][j] = dux[i][j] + duy[i][j];
  


   return res;

}


int WEN05::reconstructionX(int ind, std::vector<std::vector<double>> &ui, std::vector<std::vector<double>> &fhat){

  int i, j;
  double umm,um,u,up,upp;
  double d0n,d1n,d2n,B0n,B1n,B2n,alpha0n,alpha1n,alpha2n,alphasumn,w0n,w1n,w2n,p0n,p1n,p2n;
  d0n = 1.0/10.0; 
  d1n = 6.0/10.0; 
  d2n = 3.0/10.0; 
  double epsilon = 1e-6;
  fhat = ui;

// # pragma omp parallel for collapse(2) private(i,j,p0n,p1n,p2n,B0n,B1n,B2n,alpha0n,alpha1n,alpha2n,alphasumn,w0n,w1n,w2n,umm,um,u,up,upp)
  for(i = 0; i < Ny; i++){
    for(j = 0; j < Nx; j++){
      // Determine location of discretization
    
      // Left Boundary
      if(j==0){ 
        // if(p2d.pcoords[
        umm = p2d.vrecvWpp[i];
        um = p2d.vrecvWp[i];
        u = ui[i][j];
        up  = ui[i][j+1];
        upp = ui[i][j+2];
      }
      else if(j==1){
        umm = p2d.vrecvWp[i];
        um  = ui[i][j-1];
        u = ui[i][j];
        up  = ui[i][j+1];
        upp = ui[i][j+2];
      }
      // Right Boundary 
      else if(j==Nx-2){
        umm = ui[i][j-2];
        um  = ui[i][j-1];
        u = ui[i][j];
        up  = ui[i][j+1];
        upp = p2d.vrecvEp[i]; 
      }
      else if(j==Nx-1){
        umm = ui[i][j-2];
        um  = ui[i][j-1];
        u = ui[i][j];
        up  = p2d.vrecvEp[i];
        upp = p2d.vrecvEpp[i];
        
      } 
      // Inner Domain        
      else{
        umm = ui[i][j-2*ind];
        um  = ui[i][j-ind];
        u = ui[i][j];
        up  = ui[i][j+ind];
        upp = ui[i][j+2*ind];
      }

      // Call Reconstruction Method (Defines polynomial, beta factors etc..)
      p0n = (2.0*umm - 7.0*um + 11.0*u)/6.0;
      p1n = ( -1.0*um  + 5.0*u  + 2.0*up)/6.0;
      p2n = (2.0*u   + 5.0*up - 1.0*upp )/6.0;
      // Smooth Indicators (Beta factors)
      B0n = 13.0/12.0*pow(umm-2.0*um +u ,2.0)  + 1.0/4.0*pow(umm-4.0*um+3.0*u,2.0); 
      B1n = 13.0/12.0*pow(um -2.0*u  +up ,2.0) + 1.0/4.0*pow(um-1.0*up,2.0);
      B2n = 13.0/12.0*pow(u  -2.0*up +upp,2.0) + 1.0/4.0*pow(3.0*u-4.0*up+upp,2.0);

      // Alpha weights 
      alpha0n = d0n/pow(epsilon + B0n,2.0);
      alpha1n = d1n/pow(epsilon + B1n,2.0);
      alpha2n = d2n/pow(epsilon + B2n,2.0);

      alphasumn = alpha0n + alpha1n + alpha2n;

      // ENO stencils weigths
      w0n = alpha0n/alphasumn;
      w1n = alpha1n/alphasumn;  
      w2n = alpha2n/alphasumn;
      // Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
      fhat[i][j] = w0n*p0n + w1n*p1n + w2n*p2n;
    }
  } 

  return 0;

}

int WEN05::reconstructionY(int ind, std::vector<std::vector<double>> &ui, std::vector<std::vector<double>> &fhat){

  fhat = ui;
  int i,j;
  double umm,um,u,up,upp;
  double d0n,d1n,d2n,B0n,B1n,B2n,alpha0n,alpha1n,alpha2n,alphasumn,w0n,w1n,w2n,p0n,p1n,p2n;
  d0n = 1.0/10.0; 
  d1n = 6.0/10.0; 
  d2n = 3.0/10.0; 
  double epsilon = 1e-6;

  // # pragma omp parallel for collapse(2) private(i,j,p0n,p1n,p2n,B0n,B1n,B2n,alpha0n,alpha1n,alpha2n,alphasumn,w0n,w1n,w2n,umm,um,u,up,upp)
  for(i = 0; i < Ny; i++){
    for(j = 0; j < Nx; j++){

      // Determine location of discretization
      
      // Top Boundary 
      if(i==0){
        umm = p2d.vrecvNpp[j];
        um  = p2d.vrecvNp[j];
        u = ui[i][j];
        up  = ui[i+1][j];
        upp = ui[i+2][j];
      }
      else if(i==1){
        umm = p2d.vrecvNp[j];
        um  = ui[i-1][j];
        u = ui[i][j];
        up  = ui[i+1][j];
        upp = ui[i+2][j];
      }
      // Bottom Boundary 
      else if(i==Ny-2){
        umm = ui[i-2][j];
        um  = ui[i-1][j];
        u = ui[i][j];
        up  = ui[i+1][j];
        upp = p2d.vrecvSp[j];
      }
      else if(i==Ny-1){
        umm = ui[i-2][j];
        um  = ui[i-1][j];
        u = ui[i][j];
        up  = p2d.vrecvSp[j];
        upp = p2d.vrecvSpp[j];       
      }
        // Inner Domain        
      else{
        umm = ui[i-2*ind][j];
        um  = ui[i-ind][j];
        u = ui[i][j];
        up  = ui[i+ind][j];
        upp = ui[i+2*ind][j];
      }
      
      // Call Reconstruction Method (Defines polynomial, beta factors etc..)
      p0n = (2.0*umm - 7.0*um + 11.0*u)/6.0;
      p1n = ( -1.0*um  + 5.0*u  + 2.0*up)/6.0;
      p2n = (2.0*u   + 5.0*up - 1.0*upp )/6.0;

      // Smooth Indicators (Beta factors)
      B0n = 13.0/12.0*pow(umm-2.0*um +u ,2.0)  + 1.0/4.0*pow(umm-4.0*um+3.0*u,2.0); 
      B1n = 13.0/12.0*pow(um -2.0*u  +up ,2.0) + 1.0/4.0*pow(um-1.0*up,2.0);
      B2n = 13.0/12.0*pow(u  -2.0*up +upp,2.0) + 1.0/4.0*pow(3.0*u-4.0*up+upp,2.0);

      // Alpha weights 
      alpha0n = d0n/pow(epsilon + B0n,2.0);
      alpha1n = d1n/pow(epsilon + B1n,2.0);
      alpha2n = d2n/pow(epsilon + B2n,2.0);

      alphasumn = alpha0n + alpha1n + alpha2n;

      // ENO stencils weigths
      w0n = alpha0n/alphasumn;
      w1n = alpha1n/alphasumn;  
      w2n = alpha2n/alphasumn;

      // Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
      fhat[i][j] = w0n*p0n + w1n*p1n + w2n*p2n;
    }
  }
  return 0;
}

