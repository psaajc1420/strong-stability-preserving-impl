function unew = SSP4_step(fcn, y, h)
% usage: unew = SSP4_step(fcn, y, h)
%
% Strong Stability Preserving method of order 4 for one step of the 
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

% Inital Step
u = y;

% First Stage
f1  = feval(fcn,0,u);
u1  = u + 0.391752226571890*h*f1;

% Second Stage
f2  = feval(fcn,0,u1);
u2  = 0.444370493651235*u + 0.555629506348765*u1 + ...
      + 0.368410593050371*h*f2;

% Third Stage
f3  = feval(fcn,0,u2);
u3  = 0.620101851488403*u + 0.379898148511597*u2 + ...
    + 0.251891774271694*h*f3;

% Fourth Stage
f4  = feval(fcn,0,u3);
u4  = 0.178079954393132*u +  0.821920045606868*u3 + ... 
     + 0.544974750228521*h*f4;

% Fifth Stage
f5   =  feval(fcn,0,u4);
unew =  0.517231671970585*u2 +  0.096059710526147*u3 + ... 
      + 0.063692468666290*h*f4 + 0.386708617503269*u4 + ...
      + 0.226007483236906*h*f5;
end