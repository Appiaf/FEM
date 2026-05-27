clc;close all;
par.a = 0;  par.b = 1; par.eps = 0.5;
Nvalue = [9,99,999,9999];
par.coeffOpt = 1;
par.testOption = 1;  
par.debug = 0;
%%
l2norm   = zeros(length(Nvalue),1);
linfnorm = zeros(length(Nvalue),1);
enorm    = zeros(length(Nvalue),1);

xq  = cell(length(Nvalue),1);
error  = cell(length(Nvalue),1);
Uexact  = cell(length(Nvalue),1);
UI  = cell(length(Nvalue),1);
for j = 1:length(Nvalue)
    par.N = Nvalue(j);
    gpar = twoPointBVP(par);
    xq{j}    = gpar.xq;
    error{j} = gpar.eh_vec;
    Uexact{j} = gpar.u_vec;
    UI{j} = gpar.uh_vec;
    l2norm(j)   = gpar.L2Norm;
    linfnorm(j)= gpar.LinfNorm;
    enorm(j)    = gpar.ENorm;
end
%% make plots and convergence study
if par.eps == 0.5
figure;
plot(xq{1}, Uexact{1},'r-*','LineWidth', 1)
hold on
plot(xq{1}, UI{1},'b-o','LineWidth', 1)
hold off
grid on;
xlim([0 1]);
legend('Exact solution', 'FEM solution', 'Location','best')
xlabel('$x$', 'Interpreter', 'latex');
title(sprintf('Exact solution vs Approximate solution: $N+1=%d$, $\\epsilon=%g$', Nvalue(1)+1,par.eps), ...
      'Interpreter','latex')
filename = sprintf('img/solplot_N%d_eps%g.eps',Nvalue(1)+1,par.eps);
print('-depsc2', filename);
end

if par.eps == 0.01
figure;
plot(xq{1}, Uexact{1},'r-*','LineWidth', 1)
hold on
plot(xq{1}, UI{1},'b-o','LineWidth', 1)
hold off
grid on;
xlim([0 1]);
legend('Exact solution', 'FEM solution', 'Location','best')
xlabel('$x$', 'Interpreter', 'latex');
%ylabel('$|e_h(x)|$', 'Interpreter', 'latex');
title(sprintf('Exact solution vs Approximate solution: $N+1=%d$, $\\epsilon=%g$', Nvalue(1)+1,par.eps), ...
      'Interpreter','latex')
filename = sprintf('img/solplot_N%d_eps%g.eps',Nvalue(1)+1,par.eps);
print('-depsc2', filename);

figure;
plot(xq{2}, Uexact{2},'r-*','LineWidth', 1)
hold on
plot(xq{2}, UI{2},'b-o','LineWidth', 1)
hold off
grid on;
xlim([0 1]);
legend('Exact solution', 'FEM solution', 'Location','best')
xlabel('$x$', 'Interpreter', 'latex');
title(sprintf('Exact solution vs Approximate solution: $N+1=%d$, $\\epsilon=%g$', Nvalue(2)+1,par.eps), ...
      'Interpreter','latex')
filename = sprintf('img/solplot_N%d_eps%g.eps',Nvalue(2)+1,par.eps);
print('-depsc2', filename);
end
%% compute convergence
l2conv = zeros(length(Nvalue)-1,1);
linfconv = zeros(length(Nvalue)-1,1);
Econv = zeros(length(Nvalue)-1,1);
for j = 1 : length(Nvalue)-1
    l2conv(j) = log(l2norm(j)/l2norm(j+1))/log(10);
    linfconv(j) = log(linfnorm(j)/linfnorm(j+1))/log(10);
    Econv(j) = log(enorm(j)/enorm(j+1))/log(10);
end
