%% Strong Stability Preserving Methods with WENO

%% Fifth Order Weno Convergence Plot

% remove all figures and clear variables and screen
clc
clear
close all
figure('Units','characters','Position',[120 120 120 120]) 

nx = load('nxvals.txt');
err = load('error.txt');
loglog(nx,err,'-c*',nx,200*nx.^-5,'sy-')
title('\textbf{Convergence}','interpreter','latex','FontSize',18);
xlabel('\textbf{N}','interpreter','latex','FontSize',14);
ylabel({'\textbf{Error}'},'interpreter','latex','FontSize',14);
l = legend('\textbf{5th Order WENO}','\textbf{5th Order True}');
set(l,'Interpreter','Latex','FontSize',12)

%% SSP Convergence Plot

% clear variables and screen
clc
clear

figure('Units','characters','Position',[120 120 120 120]) 

% SSP2(2,2)
nt = load('ntvals_ssp2.txt');
err = load('error_ssp2.txt');
loglog(nt,err,'-r*',nt,0.5*nt.^-2,'g-s')
title('\textbf{Convergence}','interpreter','latex','FontSize',18);
xlabel('\textbf{N}','interpreter','latex','FontSize',14);
ylabel({'\textbf{Error}'},'interpreter','latex','FontSize',14);


hold on

% SSP3(3,3)
nt = load('ntvals_ssp3.txt');
err = load('error_ssp3.txt');
loglog(nt,err,'-bo',nt,0.1*nt.^-3,'k-d')
title('\textbf{Convergence}','interpreter','latex','FontSize',18);
xlabel('\textbf{N}','interpreter','latex','FontSize',14);
ylabel({'\textbf{Error}'},'interpreter','latex','FontSize',14);


hold on

% SSP4(4,5)
nt = load('ntvals_ssp4.txt');
err = load('error_ssp4.txt');
loglog(nt,err,'-cx',nt,0.01*nt.^-4,'y-p')
title('\textbf{Convergence}','interpreter','latex','FontSize',18);
xlabel('\textbf{N}','interpreter','latex','FontSize',14);
ylabel({'\textbf{Error}'},'interpreter','latex','FontSize',14);
l = legend('\textbf{SSP4(2,2)}','\textbf{True-2nd Order}','\textbf{SSP4(3,3)}'...
    ,'\textbf{True-3rd Order}','\textbf{SSP4(4,5)}','\textbf{True-4th Order}');
set(l,'Interpreter','Latex','FontSize',12)


%% SSP Stability Linear Videos & Final Plots

% clear variables and screen
clc
clear

% load data from txt files
x = load('x_Nonlinear.txt');
u0 = load('y0_Nonlinear.txt');

u_ssp2_CFL_09 = load('ycalc_NonlinearL_ssp2_CFL(1).txt');
u_ssp2_CFL_05 = load('ycalc_NonlinearL_ssp2_CFL(0.5).txt');
u_ssp2_CFL_15 = load('ycalc_NonlinearL_ssp2_CFL(1.5).txt');
u_erk2_CFL_09 = load('ycalc_NonlinearL_erk2_CFL(1).txt');
u_erk2_CFL_05 = load('ycalc_NonlinearL_erk2_CFL(0.5).txt');
u_erk2_CFL_15 = load('ycalc_NonlinearL_erk2_CFL(1.5).txt');

u_ssp3_CFL_09 = load('ycalc_NonlinearL_ssp3_CFL(1).txt');
u_ssp3_CFL_05 = load('ycalc_NonlinearL_ssp3_CFL(0.5).txt');
u_ssp3_CFL_15 = load('ycalc_NonlinearL_ssp3_CFL(1.5).txt');
u_erk3_CFL_09 = load('ycalc_NonlinearL_erk3_CFL(1).txt');
u_erk3_CFL_05 = load('ycalc_NonlinearL_erk3_CFL(0.5).txt');
u_erk3_CFL_15 = load('ycalc_NonlinearL_erk3_CFL(1.5).txt');

u_ssp4_CFL_14 = load('ycalc_NonlinearL_ssp4_CFL(1.4).txt');
u_ssp4_CFL_10 = load('ycalc_NonlinearL_ssp4_CFL(1).txt');
u_ssp4_CFL_20 = load('ycalc_NonlinearL_ssp4_CFL(2).txt');
u_erk4_CFL_14 = load('ycalc_NonlinearL_erk4_CFL(1.4).txt');
u_erk4_CFL_10 = load('ycalc_NonlinearL_erk4_CFL(1).txt');
u_erk4_CFL_20 = load('ycalc_NonlinearL_erk4_CFL(2).txt');


% Final Plots
figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_09(length(u_ssp2_CFL_09(:,1)),:),'-r*',x,u_erk2_CFL_09(length(u_ssp2_CFL_09(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_15(length(u_ssp2_CFL_15(:,1)),:),'-r*',x,u_erk2_CFL_15(length(u_ssp2_CFL_15(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_05(length(u_ssp2_CFL_05(:,1)),:),'-r*',x,u_erk2_CFL_05(length(u_ssp2_CFL_05(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_09(length(u_ssp3_CFL_09(:,1)),:),'-r*',x,u_erk3_CFL_09(length(u_ssp3_CFL_09(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_15(length(u_ssp3_CFL_15(:,1)),:),'-r*',x,u_erk3_CFL_15(length(u_ssp3_CFL_15(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_05(length(u_ssp3_CFL_05(:,1)),:),'-r*',x,u_erk3_CFL_05(length(u_ssp3_CFL_05(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_14(length(u_ssp4_CFL_14(:,1)),:),'-r*',x,u_erk4_CFL_14(length(u_ssp4_CFL_14(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_20(length(u_ssp4_CFL_20(:,1)),:),'-r*',x,u_erk4_CFL_20(length(u_ssp4_CFL_20(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_10(length(u_ssp4_CFL_10(:,1)),:),'-r*',x,u_erk4_CFL_10(length(u_ssp4_CFL_10(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)



%%

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_09(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_09(i,:),'-r*',x,u_erk2_CFL_09(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_15(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_15(i,:),'-r*',x,u_erk2_CFL_15(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_05(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_05(i,:),'-r*',x,u_erk2_CFL_05(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_09(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_09(i,:),'-r*',x,u_erk3_CFL_09(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_05(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_05(i,:),'-r*',x,u_erk3_CFL_05(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_15(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_15(i,:),'-r*',x,u_erk3_CFL_15(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end



figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_14(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_14(i,:),'-r*',x,u_erk4_CFL_14(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_10(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_10(i,:),'-r*',x,u_erk4_CFL_10(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end



figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_20(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_20(i,:),'-r*',x,u_erk4_CFL_20(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end



%% SSP Stability Burgers Videos & Final Plots

% clear variables and screen
clc
clear

% load data from txt files
x = load('x_Nonlinear.txt');
u0 = load('y0_Nonlinear.txt');

u_ssp2_CFL_09 = load('ycalc_NonlinearBR_ssp2_CFL(1).txt');
u_ssp2_CFL_05 = load('ycalc_NonlinearBR_ssp2_CFL(0.5).txt');
u_ssp2_CFL_15 = load('ycalc_NonlinearBR_ssp2_CFL(1.5).txt');
u_erk2_CFL_09 = load('ycalc_NonlinearBR_erk2_CFL(1).txt');
u_erk2_CFL_05 = load('ycalc_NonlinearBR_erk2_CFL(0.5).txt');
u_erk2_CFL_15 = load('ycalc_NonlinearBR_erk2_CFL(1.5).txt');

u_ssp3_CFL_09 = load('ycalc_NonlinearBR_ssp3_CFL(1).txt');
u_ssp3_CFL_05 = load('ycalc_NonlinearBR_ssp3_CFL(0.5).txt');
u_ssp3_CFL_15 = load('ycalc_NonlinearBR_ssp3_CFL(1.5).txt');
u_erk3_CFL_09 = load('ycalc_NonlinearBR_erk3_CFL(1).txt');
u_erk3_CFL_05 = load('ycalc_NonlinearBR_erk3_CFL(0.5).txt');
u_erk3_CFL_15 = load('ycalc_NonlinearBR_erk3_CFL(1.5).txt');

u_ssp4_CFL_14 = load('ycalc_NonlinearBR_ssp4_CFL(1.4).txt');
u_ssp4_CFL_10 = load('ycalc_NonlinearBR_ssp4_CFL(1).txt');
u_ssp4_CFL_20 = load('ycalc_NonlinearBR_ssp4_CFL(2).txt');
u_erk4_CFL_14 = load('ycalc_NonlinearBR_erk4_CFL(1.4).txt');
u_erk4_CFL_10 = load('ycalc_NonlinearBR_erk4_CFL(1).txt');
u_erk4_CFL_20 = load('ycalc_NonlinearBR_erk4_CFL(2).txt');


% Final Plots
figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_09(length(u_ssp2_CFL_09(:,1)),:),'-r*',x,u_erk2_CFL_09(length(u_ssp2_CFL_09(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_15(length(u_ssp2_CFL_15(:,1)),:),'-r*',x,u_erk2_CFL_15(length(u_ssp2_CFL_15(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_05(length(u_ssp2_CFL_05(:,1)),:),'-r*',x,u_erk2_CFL_05(length(u_ssp2_CFL_05(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_09(length(u_ssp3_CFL_09(:,1)),:),'-r*',x,u_erk3_CFL_09(length(u_ssp3_CFL_09(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_15(length(u_ssp3_CFL_15(:,1)),:),'-r*',x,u_erk3_CFL_15(length(u_ssp3_CFL_15(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_05(length(u_ssp3_CFL_05(:,1)),:),'-r*',x,u_erk3_CFL_05(length(u_ssp3_CFL_05(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_14(length(u_ssp4_CFL_14(:,1)),:),'-r*',x,u_erk4_CFL_14(length(u_ssp4_CFL_14(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_20(length(u_ssp4_CFL_20(:,1)),:),'-r*',x,u_erk4_CFL_20(length(u_ssp4_CFL_20(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_10(length(u_ssp4_CFL_10(:,1)),:),'-r*',x,u_erk4_CFL_10(length(u_ssp4_CFL_10(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)



%%

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_09(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_09(i,:),'-r*',x,u_erk2_CFL_09(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_15(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_15(i,:),'-r*',x,u_erk2_CFL_15(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_05(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_05(i,:),'-r*',x,u_erk2_CFL_05(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_09(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_09(i,:),'-r*',x,u_erk3_CFL_09(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_05(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_05(i,:),'-r*',x,u_erk3_CFL_05(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_15(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_15(i,:),'-r*',x,u_erk3_CFL_15(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end



figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_14(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_14(i,:),'-r*',x,u_erk4_CFL_14(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_10(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_10(i,:),'-r*',x,u_erk4_CFL_10(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end



figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_20(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_20(i,:),'-r*',x,u_erk4_CFL_20(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


%% SSP Stability Buckey-Leverett Videos & Final Plots

% clear variables and screen
clc
clear

% load data from txt files
x = load('x_Nonlinear.txt');
u0 = load('y0_Nonlinear.txt');

u_ssp2_CFL_09 = load('ycalc_NonlinearBL_ssp2_CFL(1).txt');
u_ssp2_CFL_05 = load('ycalc_NonlinearBL_ssp2_CFL(0.5).txt');
u_ssp2_CFL_15 = load('ycalc_NonlinearBL_ssp2_CFL(1.5).txt');
u_erk2_CFL_09 = load('ycalc_NonlinearBL_erk2_CFL(1).txt');
u_erk2_CFL_05 = load('ycalc_NonlinearBL_erk2_CFL(0.5).txt');
u_erk2_CFL_15 = load('ycalc_NonlinearBL_erk2_CFL(1.5).txt');

u_ssp3_CFL_09 = load('ycalc_NonlinearBL_ssp3_CFL(1).txt');
u_ssp3_CFL_05 = load('ycalc_NonlinearBL_ssp3_CFL(0.5).txt');
u_ssp3_CFL_15 = load('ycalc_NonlinearBL_ssp3_CFL(1.5).txt');
u_erk3_CFL_09 = load('ycalc_NonlinearBL_erk3_CFL(1).txt');
u_erk3_CFL_05 = load('ycalc_NonlinearBL_erk3_CFL(0.5).txt');
u_erk3_CFL_15 = load('ycalc_NonlinearBL_erk3_CFL(1.5).txt');

u_ssp4_CFL_14 = load('ycalc_NonlinearBL_ssp4_CFL(1.4).txt');
u_ssp4_CFL_10 = load('ycalc_NonlinearBL_ssp4_CFL(1).txt');
u_ssp4_CFL_20 = load('ycalc_NonlinearBL_ssp4_CFL(2).txt');
u_erk4_CFL_14 = load('ycalc_NonlinearBL_erk4_CFL(1.4).txt');
u_erk4_CFL_10 = load('ycalc_NonlinearBL_erk4_CFL(1).txt');
u_erk4_CFL_20 = load('ycalc_NonlinearBL_erk4_CFL(2).txt');


% Final Plots
figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_09(length(u_ssp2_CFL_09(:,1)),:),'-r*',x,u_erk2_CFL_09(length(u_ssp2_CFL_09(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_15(length(u_ssp2_CFL_15(:,1)),:),'-r*',x,u_erk2_CFL_15(length(u_ssp2_CFL_15(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp2_CFL_05(length(u_ssp2_CFL_05(:,1)),:),'-r*',x,u_erk2_CFL_05(length(u_ssp2_CFL_05(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_09(length(u_ssp3_CFL_09(:,1)),:),'-r*',x,u_erk3_CFL_09(length(u_ssp3_CFL_09(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_15(length(u_ssp3_CFL_15(:,1)),:),'-r*',x,u_erk3_CFL_15(length(u_ssp3_CFL_15(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp3_CFL_05(length(u_ssp3_CFL_05(:,1)),:),'-r*',x,u_erk3_CFL_05(length(u_ssp3_CFL_05(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_14(length(u_ssp4_CFL_14(:,1)),:),'-r*',x,u_erk4_CFL_14(length(u_ssp4_CFL_14(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_20(length(u_ssp4_CFL_20(:,1)),:),'-r*',x,u_erk4_CFL_20(length(u_ssp4_CFL_20(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)

figure('Units','characters','Position',[0 0 120 120]) 
plot(x,u0,'-',x,u_ssp4_CFL_10(length(u_ssp4_CFL_10(:,1)),:),'-r*',x,u_erk4_CFL_10(length(u_ssp4_CFL_10(:,1)),:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)


%%

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_09(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_09(i,:),'-r*',x,u_erk2_CFL_09(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_15(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_15(i,:),'-r*',x,u_erk2_CFL_15(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp2_CFL_05(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp2_CFL_05(i,:),'-r*',x,u_erk2_CFL_05(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(2,2)','ERK2');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_09(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_09(i,:),'-r*',x,u_erk3_CFL_09(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_05(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_05(i,:),'-r*',x,u_erk3_CFL_05(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end

figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp3_CFL_15(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp3_CFL_15(i,:),'-r*',x,u_erk3_CFL_15(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(3,3)','ERK3');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end



figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_14(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_14(i,:),'-r*',x,u_erk4_CFL_14(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end


figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_10(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_10(i,:),'-r*',x,u_erk4_CFL_10(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end



figure('Units','characters','Position',[0 0 120 120]) 
for i = 1:length(u_ssp4_CFL_20(:,1))
% Plot Solutions
plot(x,u0,'-',x,u_ssp4_CFL_20(i,:),'-r*',x,u_erk4_CFL_20(i,:),'-bs');
axis([-1 1 .9 2.1])
title('\textbf{Stability of Methods}','interpreter','latex','FontSize',18);
xlabel('\textbf{x}','interpreter','latex','FontSize',14);
ylabel({'\textbf{u(x)}'},'interpreter','latex','FontSize',14);
l = legend('Initial Condtion','SSP(4,5)','ERK4');
set(l,'Interpreter','Latex','FontSize',12)
pause(0.2);
end





