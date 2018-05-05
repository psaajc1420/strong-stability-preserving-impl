function unew = RK4_step(fcn, t, u, h)
% usage: unew = RK4_step(fcn, t, u, h)
%
% Runge-Kutta method of order 4 for one step of the scalar-valued
% ODE problem, 
%     y' = f(t,y), t in tspan,
%     y(t0) = y0.
%
% Inputs:  fcn = function name for ODE right-hand side, f(t,y)
%          t = current time
%          u = current solution
%          h = time step size
% Outputs: unew = updated solution
%
% Jacob Cadena
% Math 6321, SMU
% Fall 2016

% set coefficients
b = [1, 2, 2, 1]/6;
c = [0, 1/2, 1/2, 1];
A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];

% compute z1 and k1
z1 = u;
k1 = feval(fcn, t+c(1)*h, z1);

% compute z2 and k2
z2 = u + h*A(2,1)*k1;
k2 = feval(fcn, t+c(2)*h, z2);

% compute z3 and k3
z3 = u + h*(A(3,1)*k1 + A(3,2)*k2);
k3 = feval(fcn, t+c(3)*h, z3);

z4 = u + h*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);
k4 = feval(fcn, t+c(4)*h, z4);
% update solution in time
unew = u + h*(b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4);

end