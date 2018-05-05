function unew = SSP3_step(fcn, y, h)
% usage: unew = SSP3_step(fcn, y, h)
%
% Strong Stability Preserving method of order 3 for one step of the 
% scalar-valued ODE problem, 
%     y' = f(y), 
%     y(t0) = y0.
%
% Inputs:  fcn = function name for ODE right-hand side, f(y)
%          y = current solution
%          h = time step size
% Outputs: unew = updated solution
%
% Jacob Cadena
% Math 6321, SMU
% Fall 2016

u = y;

f1 = feval(fcn,0,u);
u1 = u + h*f1;

f2 = feval(fcn,0,u1);
u2 = (3*u + (u1 + h *f2))/4;

f3 = feval(fcn,0,u2);
unew = (u + 2*u2 + 2*h*f3)/3;

end