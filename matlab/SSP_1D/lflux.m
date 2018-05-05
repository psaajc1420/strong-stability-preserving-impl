function fhat = lflux(u,ind)

% Shifted for stencil
umm = circshift(u,2*ind);
um  = circshift(u,ind);
up  = circshift(u,-ind);
upp = circshift(u,-2*ind);

% Polynomials
p0p = ( -umm + 5*um + 2*u  )/6;
p1p = ( 2*um + 5*u  - up   )/6;
p2p = (11*u  - 7*up + 2*upp)/6;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

for i = 1:length(u)
% Smooth Indicators (Beta factors)
B0p = 13/12*(umm(i)-2*um(i)+u(i)  )^2 + 1/4*(umm(i)-4*um(i)+3*u(i))^2; 
B1p = 13/12*(um(i) -2*u(i) +up(i) )^2 + 1/4*(um(i)-up(i))^2;
B2p = 13/12*(u(i)  -2*up(i)+upp(i))^2 + 1/4*(3*u(i) -4*up(i)+upp(i))^2;

% Alpha weights 
alpha0p = d0p/(epsilon + B0p)^2;
alpha1p = d1p/(epsilon + B1p)^2;
alpha2p = d2p/(epsilon + B2p)^2;

alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p/alphasump;
w1p = alpha1p/alphasump;
w2p = alpha2p/alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
fhat(i) = w0p*p0p(i) + w1p*p1p(i) + w2p*p2p(i);

end

end