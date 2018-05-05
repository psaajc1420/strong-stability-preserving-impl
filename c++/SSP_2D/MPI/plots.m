%% Plotting results

%% clear screen
clc
clear

%% Weak Scaling

t = [2.27 2.5 2.63 3.36];
N = [100 200 400 800];

figure()
plot(N,t);
axis([100 800 0 5])
title('$\bf{Weak \;Scaling}$','Interpreter','Latex','FontSize',18);
xlabel('$\textbf{P}$','Interpreter','Latex','FontSize',14)
ylabel('$\textbf{T}$','Interpreter','Latex','FontSize',14)



%% Strong Scaling

t = [61.58 32.31 18 8.52 4.27 2.23 1.17 .71 .42];
p = [1 2 4 8 16 32 64 128 256];

figure()
loglog(p,t);
title('$\bf{Strong \;Scaling}$','Interpreter','Latex','FontSize',18);
xlabel('$\bf{\log(p)}$','Interpreter','Latex','FontSize',14)
ylabel('$\bf{\log(T)}$','Interpreter','Latex','FontSize',14)


%% Write to text file

fileID = fopen('time_data.txt','w');
fprintf(fileID,'%.2f\n',t);
fclose(fileID);
