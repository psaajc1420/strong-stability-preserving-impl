%% SSP Methods with Burgers Flux

%% Clear and Close All Windows
clear; close all; clc;

%% Delete Files
delete *.avi 

%% Parameters

% Build discrete domain
nx = 160;
a=-1; b=1; dx=(b-a)/nx; x=a+dx:dx:b;   
T = 1;                                 % Final Time                                 
t0 = 0;                                % Start Time

% Linear Flux
flux = @(w) w; 
dflux = @(w) ones(size(w)); 
fprintf('Linear Flux:\n'); 

% Source term
S = @(w) zeros(size(w));

% Intial Condition
xmid=0.5*(x(end)+x(1));
u0 = ones(size(x));
u0(x<=xmid) = 2;


f = @(t,u) (WENO5(u,flux,dflux,S,dx));


%% Order 2 Method Stability Preserving Demo
CFLS = [1 1.5 0.5];
for j = 1:length(CFLS)
    
U_ssp2 = u0; U_rk2 = u0;
dt = CFLS(j)*dx/max(abs(U_ssp2));
nt = (T-t0)/dt;
tvals = linspace(t0,T,nt);
figure('Units','characters','Position',[0 0 120 120])        
writer = VideoWriter(sprintf('SSP2-RK2_CFL(%g).avi',CFLS(j)));
open(writer);

    for i = 1:nt-1
    % Update/correct time step
    t = tvals(i);
    % Solve odes
    U_ssp2 = SSP2_step(f, U_ssp2, dt);
    U_rk2 = RK2_step(f, t, U_rk2, dt);

    % Plot Solutions
    plot(x,u0,'-',x,U_ssp2,'-r*',x,U_rk2,'-bs');
    axis([-1 1 .8 2.2])
    title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
    xlabel('\textbf{x}','interpreter','latex','FontSize',14);
    ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
    l = legend('Initial Condition','SSP(2,2)','RK2');
    set(l,'Interpreter','Latex','FontSize',12);

    % Store the frame
    frame = getframe(gcf); 
    writeVideo(writer,frame);

    end        

close(writer);

end
    
%% Order 3 Method Stability Preserving Demo
CFLS = [1 2.5 0.5];
for j = 1:length(CFLS)
U_ssp3 = u0; U_rk3 = u0; 
dt = CFLS(j)*dx/max(abs(U_ssp3));
nt = (T-t0)/dt;
tvals = linspace(t0,T,nt);
figure('Units','characters','Position',[0 0 120 120])  
writer = VideoWriter(sprintf('SSP3-RK3_CFL(%g).avi',CFLS(j)));
open(writer);

    for i = 1:nt-1
    % Update time
    t = tvals(i);
    % Solve odes
    U_ssp3 = SSP3_step(f, U_ssp3, dt);
    U_rk3 = RK3_step(f, t, U_rk3, dt);

    % Plot Solutions
    plot(x,u0,'-',x,U_ssp3,'-r*',x,U_rk3,'-bs');
    axis([-1 1 .9 2.1])
    title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
    xlabel('\textbf{x}','interpreter','latex','FontSize',14);
    ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
    l = legend('Initial Condtion','SSP(3,3)','RK3');
    set(l,'Interpreter','Latex','FontSize',12)

    % Store the frame
    frame = getframe(gcf); 
    writeVideo(writer,frame);

    end

close(writer);
end

%% Order 4 Method Stability Preserving Demo
CFLS = [1.4 2 0.8];
for j = 1:length(CFLS)   
U_ssp4 = u0; U_rk4 = u0;
dt = CFLS(j)*dx/max(abs(U_ssp4));
nt = (T-t0)/dt;
tvals = linspace(t0,T,nt);
figure('Units','characters','Position',[0 0 120 120])  
writer = VideoWriter(sprintf('SSP4-RK4_CFL(%g).avi',CFLS(j)));
open(writer);

    for i = 1:nt-1
    % Update time
    t = tvals(i);
    % Solve odes
    U_ssp4 = SSP4_step(f, U_ssp4, dt);
    U_rk4 = RK4_step(f, t, U_rk4, dt);

    % Plot Solutions
    plot(x,u0,'-',x,U_ssp4,'-r*',x,U_rk4,'-bs');
    axis([-1 1 .9 2.1])
    title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
    xlabel('\textbf{x}','interpreter','latex','FontSize',14);
    ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
    l = legend('Initial Condition','SSP(5,4)','RK4');
    set(l,'Interpreter','Latex','FontSize',12)

    % Store the frame
    frame = getframe(gcf); 
    writeVideo(writer,frame);
    
    end

close(writer);

end

