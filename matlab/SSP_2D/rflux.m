function fhat = rflux(u,ind)

% Shifted for stencil
umm = circshift(u,ind*2);
um  = circshift(u,ind);
up  = circshift(u,-ind);
upp = circshift(u,-ind*2);

% Polynomials
p0n = (2*umm - 7*um + 11*u)/6.;
p1n = ( -um  + 5*u  + 2*up)/6.;
p2n = (2*u   + 5*up - upp )/6.;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;


% Smooth Indicators (Beta factors)
B0n= 13/12*( umm-2*um+u ).^2 + 1/4*(umm-4*um+3*u).^2; 
B1n = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
B2n = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u-4*up+upp).^2;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;

alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
fhat = w0n.*p0n + w1n.*p1n + w2n.*p2n;

end