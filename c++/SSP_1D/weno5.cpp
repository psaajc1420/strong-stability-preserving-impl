/* Fifth order WENO method file.

   Jacob Cadena
   Math 6321 @ SMU
   Fall 2016  */

// Inclusions    
#include "weno5.hpp"

// constructor (sets RHS function pointer, allocates local data)
WEN05::WEN05(std::vector<double>& u_, RHSFlux& flux_, RHSSource& source_,double dx_) { 

  flux = &flux_;
  source = &source_;
  uf = u_;
  v = u_;
  w = u_;
  s = u_;
  dx = dx_;

  N = u_.size();

  // Constants for Right FLux
  d0n = 1.0/10.0; 
  d1n = 6.0/10.0; 
  d2n = 3.0/10.0; 
  epsilon = 1e-6;

}

int WEN05::fluxSplitting(std::vector<double> &u, std::vector<double> &f, std::vector<double> &df, std::vector<double> &vp
                                   , std::vector<double> &vn){

    a = max(abs(df));

    for(int i =0; i < N; i++){
      vp[i] = 0.5*(f[i] + a*u[i]);
    }

    for(int i =0; i < N; i++){
      vn[i] = 0.5*(f[i] - a*u[i]);
    }

    return 0;
}

std::vector<double> WEN05::residual(std::vector<double> &res) {
    
     w  = res;
     vn = res;
     vp = res;
     f  = res;
     df = res;

    // Lax-Friedrichs Flux Splitting
    if (flux->Evaluate(w,f,df) != 0) {
     std::cerr << "Step: Error in Flux function\n";
     return res;
    }


    if(source->Evaluate(w,s) !=0) {
     std::cerr << "Step: Error in Source function\n";
     return res;
    }

    fluxSplitting(w, f, df, vp, vn); 


    reconstruction( 1, vp, rf);
    reconstruction(-1, vn, lf);

    res[0] =  -((lf[0] - lf[1] + rf[0] - rf[N-1])/dx - s[0]);
    for(int i =1; i < N-1; i++){
      res[i] =  -((lf[i] - lf[i+1] + rf[i] - rf[i-1])/dx - s[i]);
    }
    res[N-1] =  -((lf[N-1] - lf[0] + rf[N-1] - rf[N-2])/dx - s[N-1]);

   return res;

}



int WEN05::periodicBoundary(int i,int ind, std::vector<double>& ui){
     
      // Implementing periodic conditions
      if(i==0){
      umm = ui[N-2];
      um  = ui[N-1];
      u = ui[i];
      up  = ui[i+1];
      upp = ui[i+2];
    }
    else if(i==1){
      umm = ui[N-1];
      um  = ui[i-1];
      u = ui[i];
      up  = ui[i+1];
      upp = ui[i+2];
    }
    else if(i==N-2){
      umm = ui[i-2];
      um  = ui[i-1];
      u = ui[i];
      up  = ui[i+1];
      upp = ui[0]; 
    }
    else if(i==N-1){
      umm = ui[i-2];
      um  = ui[i-1];
      u = ui[i];
      up  = ui[0];
      upp = ui[1]; 
      
    }        
    else{
      umm = ui[i-2*ind];
      um  = ui[i-ind];
      u = ui[i];
      up  = ui[i+ind];
      upp = ui[i+2*ind];
    }
    
    return 0;
}
int WEN05::reconstruction(int ind, std::vector<double> &ui, std::vector<double> &fhat){

  fhat = ui;

  for(int i = 0; i < N; i++){

    // Determine location of discretization
    periodicBoundary(i, ind, ui);

    // Polynomials
    p0n = (2.0*umm - 7.0*um + 11.0*u)/6.0;
    p1n = ( -1.0*um  + 5.0*u  + 2.0*up)/6.0;
    p2n = (2.0*u   + 5.0*up - 1.0*upp )/6.0;
    // Smooth Indicators (Beta factors)
    B0n = 13.0/12.0*pow(umm-2.0*um +u ,2)  + 1.0/4.0*pow(umm-4.0*um+3.0*u,2); 
    B1n = 13.0/12.0*pow(um -2.0*u  +up ,2) + 1.0/4.0*pow(um-1.0*up,2);
    B2n = 13.0/12.0*pow(u  -2.0*up +upp,2) + 1.0/4.0*pow(3.0*u-4.0*up+upp,2);

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
    fhat[i] = w0n*p0n + w1n*p1n + w2n*p2n;
  }

   return 0;

}

