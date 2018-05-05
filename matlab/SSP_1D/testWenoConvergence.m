clc
clear

flux = @(w) w;
dflux = @(w) ones(size(w));
S = @(w) zeros(size(w));
nx  = 20:20:140;
nt = 1000;
tvals = linspace(0,1,nt);
dt = tvals(2) - tvals(1);


for j = 1:length(nx)
    
% Set up tests for this h values               
a=-1; b=1; dx=(b-a)/nx(j); %x=a+dx:dx:b;     % Spatial mesh size
% load initial conditions
x = linspace(a + dx,b,nx(j));
% dx = x(2) - x(1);
f = @(t,u) (WENO5(u,flux,dflux,S,dx));
u0 = sin(pi*x);
u=u0;

for i = 1:nt-1    
    u = RK3_step(f,tvals(i), u, dt);
    exact = sin(pi*(x-tvals(i+1))); 
    err = exact - u;
    error(i) = norm(err,'inf');
end
   
    err_WENO5(j) = max(abs(error));
end

err = err_WENO5;
loglog(nx,err,'-c',nx,200*nx.^-5,'y-')
fprintf('n = %10g,  err = %.2e\n',nx(1),err(1))
for i=2:length(nx)
   fprintf('n = %10g,  err = %.2e,  rate = %g\n',...
   nx(i),err(i),-log(err(i)/err(i-1))/log(nx(i)/nx(i-1)))
end
    




