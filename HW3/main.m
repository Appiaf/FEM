clc;
clear;
close all;

%% Parameters
Nvalue = [9,19,39,79,159];
%Nvalue = 3;
par.xa = 0;
par.xb = 1;
par.ya = 0;
par.yb = 1;
par.p = 3;
par.q = 2;

%par.meshType = "uniform";   % or "non_uniform"
par.meshType = "nonuniform";
par.testfunc = 4;           % use 4 for meaningful convergence

%% Storage
l2norm = zeros(length(Nvalue),1);
enorm  = zeros(length(Nvalue),1);
AssembleTime  = zeros(length(Nvalue),1);
SolveMatTime  = zeros(length(Nvalue),1);
%% Solve for each mesh
for j = 1:length(Nvalue)
    par.N = Nvalue(j);
    [errL2, errH1,MatAssemTime,solveMat] = twoPointBVP2D(par);
    l2norm(j) = errL2;
    enorm(j)  = errH1;
    AssembleTime(j) = MatAssemTime;
    SolveMatTime(j) = solveMat;
end

%% Compute mesh size h
h = 1./(Nvalue(:)-1);

%% Compute convergence rates
l2conv = NaN(length(Nvalue),1);
Econv  = NaN(length(Nvalue),1);

for j = 1:length(Nvalue)-1
    l2conv(j+1) = log(l2norm(j)/l2norm(j+1)) / log(2);
    Econv(j+1)  = log(enorm(j)/enorm(j+1)) / log(2);
end

%% Create results table
Results = table(Nvalue(:), h, l2norm, l2conv, enorm, Econv, ...
    AssembleTime,SolveMatTime,'VariableNames',...
    {'N','h','L2_Error','L2_Order','H1_Error','H1_Order','cpuTime_MatAss(s)','cpuTime_solvMat(s)'});

%% Display results
format short e
disp('======================================')
disp(['Mesh Type   : ', char(par.meshType)])
disp(['Test Function: ', num2str(par.testfunc)])
disp('======================================')
disp(Results)

%% Optional: save to file
writetable(Results, 'FEM_results.csv');

%% Optional: log-log plot
% figure
% loglog(h, l2norm, '-o', 'LineWidth',2)
% hold on
% loglog(h, enorm, '-s', 'LineWidth',2)
% grid on
% xlabel('h')
% ylabel('Error')
% legend('L2 Error','H1 Error','Location','southEast')
% title('FEM Convergence')