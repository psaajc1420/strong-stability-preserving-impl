function unew = SSP2_step(fcn, y, h)
% usage: unew = SSP2_step(fcn, y, h)
%
% Strong Stability Preserving method of order 2 for one step of the 
% scalar-valued ODE problem, 
%     y' = f(y), t in tspan,
%     y(t0) = y0.
%
% Inputs:  fcn = function name for ODE right-hand side, f(y)
%          t = current time
%          u = current solution
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
unew = 0.5*(u + u1 + h*f2);

end