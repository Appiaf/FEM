clc;close all;
par.a = 0;  par.b = 1;
Nvalue = [9,19,39,79,159,319];
par.coeffOpt = 1;
par.testOption = 1;  %1:u=x(x-1)(sin5x+3exp(x)),2:u=4,3:u=x-2,4:u=x^2-3
par.debug = 0;
%%
l2norm   = zeros(length(Nvalue),1);
linfnorm = zeros(length(Nvalue),1);
enorm    = zeros(length(Nvalue),1);

xq  = cell(length(Nvalue),1);
error  = cell(length(Nvalue),1);
for j = 1:length(Nvalue)
    par.N = Nvalue(j);
    gpar = twoPointBVP(par);
    xq{j}    = gpar.xq;
    error{j} = gpar.eh_vec;
    l2norm(j)   = gpar.L2Norm;
    linfnorm(j)= gpar.LinfNorm;
    enorm(j)    = gpar.ENorm;
end
%% make error plots and convergence study
figure; 
plot(xq{1}, error{1},'g-*','LineWidth', 1)
hold on
plot(xq{2}, error{2},'r-*','LineWidth', 1)
hold on
plot(xq{3}, error{3},'b-*','LineWidth', 1)
hold off
grid on;
xlim([0 1]);
legend('N = 10', 'N = 20', 'N = 40', 'Location','best')
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$e_h(x)$', 'Interpreter', 'latex');
title('error function graph' )
filename = sprintf('img/errorgraph.eps');
print('-depsc2', filename);


figure; 
semilogy(xq{1}, abs(error{1}),'g-*','LineWidth', 1)
hold on
semilogy(xq{2}, abs(error{2}),'r-*','LineWidth', 1)
hold on
semilogy(xq{3}, abs(error{3}),'b-*','LineWidth', 1)
hold off
grid on;
xlim([0 1]);
legend('N = 10', 'N = 20', 'N = 40', 'Location','best')
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$|e_h(x)|$', 'Interpreter', 'latex');
title('semi-log error function graph' )
filename = sprintf('img/semilogerrgraph.eps');
print('-depsc2', filename);
%% compute convergence
l2conv = zeros(length(Nvalue)-1,1);
linfconv = zeros(length(Nvalue)-1,1);
Econv = zeros(length(Nvalue)-1,1);
for j = 1 : length(Nvalue)-1
    l2conv(j) = log(l2norm(j)/l2norm(j+1))/log(2);
    linfconv(j) = log(linfnorm(j)/linfnorm(j+1))/log(2);
    Econv(j) = log(enorm(j)/enorm(j+1))/log(2);
end
