
clear; close all; clc;

cu =  1.0;  % Scalar velocity in x direction
cv =  1.0;  % Scalar velocity in y direction
nt = 1000;
tvals = linspace(0,1,nt);
dt = tvals(2) - tvals(1);
nxs = 80;
nys = 80;
% nxs = 40*2.^(0:4);
% nys = 40*2.^(0:4);

% Create Fluxes and Source
f = @(w) cu*w;
g = @(w) cv*w;
df = @(w) cu*ones(size(w));
dg = @(w) cv*ones(size(w));
S = @(w) zeros(size(w));

% Initial Condition
type = 'sines';

% Create Dimensions for Figure
fig = figure(1);
setWindowPosition(fig,500,1000)

% Loop through different number of grid pts
% to find convergence
for j = 1:length(nxs)

nx = nxs(j);
ny = nys(j);

a=0; b=1;
dx=(b-a)/nx; %x=a+dx:dx:b;     % Spatial mesh size
dy=(b-a)/ny; %y=a+dy:dy:b;     % Spatial mesh size
% load initial conditions
x = a+dx:dx:b;
y = a+dy:dy:b;

[xx,yy]=meshgrid(x,y);

q0 = initialConditions(xx,yy,type);

% load initial conditions
q=q0; it=0;

% Set the view of mesh plot
az = 56;
el = 8;

% set plot range
plotrange = [0, b/dy, 0, b/dy -1, 1];

L = @(t,q) residual(q,f,df,g,dg,dx,dy,S,'LF');

for kt=2:nt
    
    t = tvals(kt);

    true = exact(xx,yy,t,type);
    q = SSP4_step(L,q,dt);
     
    if(rem(kt,10)==0)
        subplot(1,2,1); mesh(q); colormap Copper; view(az,el); %axis(plotrange);
        %colorbar('location','EastOutside');
        title(['WENO5, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(tvals(kt))])
        xlabel('x points'); ylabel('y points'); zlabel('q(x,y)');
        axis(plotrange)
        subplot(1,2,2); mesh(true); colormap Copper; view(az,el); %axis(plotrange);% contourf(q); 
        title(['WENO5, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(tvals(kt))])
        xlabel('x points'); ylabel('y points');
        axis(plotrange)
        drawnow
        
    end

    err = true - q;
    error(kt,1) = norm(err,1);
    error(kt,2) = norm(err,'inf');  
    
end
  fprintf('N = %d, M = %d is done\n',nx,ny)
  err_WENO5(j,:) = max(abs(error));


end
%--------------------------------------------------------------
figure(2)
err = err_WENO5;
loglog(nxs,err(:,1),'-c',nxs,200*nxs.^-5,'y-')
figure(3)
loglog(nxs,err(:,2),'-c',nxs,200*nxs.^-5,'y-')
fprintf('n = %10g,  err = %.2e\n',nxs(1),err(:,1))
for i=2:length(nxs)
   fprintf('n = %10g,  err 1 = %.2e,  rate 1 = %g  err 1 = %.2e,  rate 1 = %g\n',...
   nxs(i),err(i),-log(err(i,1)/err(i-1,1))/log(nxs(i)/nxs(i-1)), ...
   err(i),-log(err(i,2)/err(i-1,2))/log(nxs(i)/nxs(i-1)))
end
