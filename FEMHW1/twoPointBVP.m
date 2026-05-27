function [gpar] = twoPointBVP(par)
a = par.a; b = par.b;
N = par.N; debug = par.debug;
dx = 1/(N + 1);
x = zeros(N+2,1);
h = zeros(N+1,1);
c = zeros(N+2,1);
coeffOpt = par.coeffOpt;
testOption = par.testOption;
% nonuniform mesh
for i = 1 : N+1
    if mod(i,2) == 0
        h(i) = 1.1*dx;
    else
        h(i) = 0.9*dx;
    end
end
x(1) = a;
for j = 1 : N+1
    x(j+1) = h(j) + x(j);
end
% define functions
if testOption == 1
    UTRUE   = @(x) x.*(x-1).*(sin(5*x)+3*exp(x));
    UTRUEX  = @(x) x*(x-1)*(5*cos(5*x) + 3*exp(x)) ...
        + (sin(5*x)+3*exp(x))*(2*x-1);
    UTRUEXX = @(x) x*(x-1)*(-25*sin(5*x) + 3*exp(x)) ...
        + 2*(5*cos(5*x)+3*exp(x))*(2*x-1) + 2*(sin(5*x)+3*exp(x));
elseif testOption == 2
    UTRUE   = @(x) 4;
    UTRUEX  = @(x) 0;
    UTRUEXX = @(x) 0;
elseif testOption == 3
    UTRUE   = @(x) x-2;
    UTRUEX  = @(x) 1;
    UTRUEXX = @(x) 0;
elseif testOption == 4
    UTRUE   = @(x) x.^2 - 3;
    UTRUEX  = @(x) 2*x;
    UTRUEXX = @(x) 2;
end
alpha = UTRUE(a);
beta  = UTRUE(b);
if coeffOpt == 1
    f   = @(x) -3*UTRUEXX(x) + 2*UTRUE(x);
else
    f   = @(x) -(1+x).*UTRUEXX(x) - UTRUEX(x);
end
% form matrix A
A = formMatrix(x,h,N,coeffOpt);
%form vector F
F_hat = zeros(N,1);
for i = 1:N
    F_hat(i) = f(x(i+2))*h(i+1)/6 + f(x(i+1))*(h(i)+h(i+1))/3 + f(x(i))*h(i)/6;
end
if coeffOpt == 1
    F_til = zeros(N,1);
    F_til(1) = 3*alpha/h(1) - 2*alpha*h(1)/6;
    F_til(N) = 3*beta/h(N+1) - 2*beta*h(N+1)/6;
    F = F_hat + F_til;
else
    F = F_hat;
end
% solve for c
c(1) = alpha;
c(N+2) = beta;
c(2:N+1) = A\F;

%% calculate errors norms
M = 3;
L2Norm = 0;
LinfNorm = 0;
ENorm = 0;
for j = 1: N+1
    for m = 1:M
        idx = x(j) + h(j)*m/M;
        e_h = UhFunc(idx,j,x,h,c) - UTRUE(idx);
        e_hx = UhFuncx(j,h,c) - UTRUEX(idx);
        L2Norm = L2Norm + e_h^2* h(j)/M;
        LinfNorm = max(LinfNorm ,abs(e_h));
        ENorm = ENorm + e_hx^2* h(j)/M;
    end
end
gpar.L2Norm = sqrt(L2Norm);
gpar.ENorm = sqrt(ENorm);
gpar.LinfNorm = LinfNorm;
%% error function
Mq = 20;
e_h0 = UhFunc(x(1),1,x,h,c) - alpha;
xq_vec = x(1);
eh_vec = e_h0;
for j = 1: N+1
    for m = 1:Mq
        idxq = x(j) + h(j)*m/Mq;
        e_h = (UhFunc(idxq,j,x,h,c,N) - UTRUE(idxq));
        xq_vec = [xq_vec;idxq];
        eh_vec = [eh_vec;e_h];
    end
end
gpar.xq = xq_vec;
gpar.eh_vec = eh_vec;
%%
if debug == 1
    % for testing uh function at the nodes
    uh = zeros(length(x),1);
    uexact = zeros(length(x),1);
    for j  = 1:N+1
        uh(j)  = UhFunc(x(j),j,x,h,c,N);
        uexact(j) = UTRUE(x(j));
    end
    uh(N+2)  = UhFunc(x(N+2),N+1,x,h,c,N);
    uexact(N+2) = UTRUE(x(N+2));
    figure(1)
    plot(x, uh,'b-*','LineWidth', 1.5);
    hold on
    plot(x, uexact,'r-','LineWidth', 1.5);
    hold on
    plot(x, c,'k--','LineWidth', 1.5);
    hold off;
    grid on;
end
return
end