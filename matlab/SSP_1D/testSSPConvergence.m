%% Order 2 Method Convergence
clc
clear
f   = @(t,y) (-y);
T   = 1;
T0  = 0;
u0  = 1;
nt  = 10*2.^(1:5);
fig = figure();
set(fig,'Units','characters','Position',[120 120 120 120])

for j = 1:length(nt)
    % Set up tests for this h values
    u_ssp2(1) = u0;
    h(j) = (T-T0)/nt(j);
    tvals = T0:h(j):T;
    ytrue = exp(-tvals);

    for i = 1:length(tvals)-1
        u_ssp2(i+1) = SSP2_step(f, u_ssp2(i), h(j));
    end
    % Compute error
    err = ytrue - u_ssp2; err_SSP2(j) = max(abs(err));
end
% output results
disp('Results:')

err = err_SSP2;
loglog(nt,err_SSP2,'r-',nt,nt.^-2,'g-')
fprintf('n = %10g,  err = %.2e\n',nt(1),err(1))
for i=2:length(nt)
   fprintf('n = %10g,  err = %.2e,  rate = %g\n',...
   nt(i),err(i),-log(err(i)/err(i-1))/log(nt(i)/nt(i-1)))
end

hold on
%% Order 3 Method Convergence
for j = 1:length(nt)
    % Set up tests for this h values
    u_ssp3(1) = u0;
    h(j) = (T-T0)/nt(j);
    tvals = T0:h(j):T;
    ytrue = exp(-tvals);

    for i = 1:length(tvals)-1
        u_ssp3(i+1) = SSP3_step(f, u_ssp3(i), h(j));
    end
    % Compute error
    err = ytrue - u_ssp3; err_SSP3(j) = max(abs(err));
end
% output results
disp('Results:')

err = err_SSP3;
loglog(nt,err_SSP3,'b-',nt,nt.^-3,'k-')
fprintf('n = %10g,  err = %.2e\n',nt(1),err(1))
for i=2:length(nt)
   fprintf('n = %10g,  err = %.2e,  rate = %g\n',...
   nt(i),err(i),-log(err(i)/err(i-1))/log(nt(i)/nt(i-1)))
end


hold on
%% Order 4 Method Convergence

for j = 1:length(nt)
    % Set up tests for this h values
    u_ssp4(1) = u0;  u_rk4(1) = u0; 
    h(j) = (T-T0)/nt(j);
    tvals = T0:h(j):T;
    ytrue = exp(-tvals);

    for i = 1:length(tvals)-1
        u_ssp4(i+1) = SSP4_step(f, u_ssp4(i), h(j));
    end
    % Compute error
    err = ytrue - u_ssp4; err_SSP4(j) = max(abs(err));
end
% output results
disp('Results:')

err = err_SSP4;
loglog(nt,err,'-c',nt,.1*nt.^-4,'y-')
fprintf('n = %10g,  err = %.2e\n',nt(1),err(1))
for i=2:length(nt)
   fprintf('n = %10g,  err = %.2e,  rate = %g\n',...
   nt(i),err(i),-log(err(i)/err(i-1))/log(nt(i)/nt(i-1)))
end
legends = {'SSP(2,2)','SL = 2','SSP(3,3)','SL = 3','SSP(5,4)','SL = 4'};

title('\textbf{Convergence for SSP Methods}','Interpreter','LaTex','FontSize',16);
xlabel('$\textbf{N}$','Interpreter','LaTex','FontSize',12);
ylabel('$\textbf{Error}$','Interpreter','LaTex','FontSize',12);
legend(legends,'Location','best','Interpreter','Latex','FontSize',12)