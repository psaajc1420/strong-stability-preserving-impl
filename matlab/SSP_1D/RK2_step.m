function ynew = RK2_step(fcn, t, y, h)
% usage: ynew = RK2_step(fcn, t, y, h)
%
% Runge-Kutta of order 2 solver for one step of the scalar-valued ODE
% problem, 
%     y' = f(t,y), t in tspan,
%     y(t0) = y0.
%
% Inputs:  fcn = function name for ODE right-hand side, f(t,y)
%          t = current time
%          y = current solution
%          h = time step size
% Outputs: ynew = updated solution
%
% Jacob Cadena
% Math 6321, SMU
% Fall 2016

% get ODE RHS at this time step
f = feval(fcn, t, y);

% set z1 and z2
z1 = y;
z2 = y + h/2*f;

% get f(t+h/2,z2)
f2 = feval(fcn, t+h/2, z2);

% update solution in time
ynew = y + h*f2;

end