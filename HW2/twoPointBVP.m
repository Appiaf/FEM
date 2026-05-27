function [gpar] = twoPointBVP(par)
a = par.a; b = par.b; eps = par.eps;
N = par.N; debug = par.debug;
h = 1/(N + 1);
x = zeros(N+2,1);
c = zeros(N+2,1);
x(1) = a;
for j = 1 : N+1
    x(j+1) = h + x(j);
end
UTRUE   = @(x) (1-(exp((x-1)/eps)))/(1-exp(-1/eps));
UTRUEX  = @(x) - 1/eps * (exp((x-1)/eps))/(1-exp(-1/eps));
alpha = 1;
beta  = 0;
A = formMatrix(x,h,par);
%form vector F
F = zeros(N,1);
F(1) = eps*alpha/h + alpha/2;
F(N) = eps*beta/h - beta/2;
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
        idx = x(j) + h*m/M;
        e_h = UhFunc(idx,j,x,h,c) - UTRUE(idx);
        e_hx = UhFuncx(j,h,c) - UTRUEX(idx);
        L2Norm = L2Norm + e_h^2* h/M;
        LinfNorm = max(LinfNorm ,abs(e_h));
        ENorm = ENorm + e_hx^2* h/M;
    end
end
gpar.L2Norm = sqrt(L2Norm);
gpar.ENorm = sqrt(ENorm);
gpar.LinfNorm = LinfNorm;
%% FEM solution function
Mq = 20;
xq_vec = x(1);
% store FEM and exact solutions
uh_vec = UhFunc(x(1),1,x,h,c);
u_vec  = UTRUE(x(1));
eh_vec = u_vec -uh_vec;
for j = 1:N+1
    for m = 1:Mq
        idxq = x(j) + h*m/Mq;
        % compute solutions
        uh = UhFunc(idxq,j,x,h,c);
        u  = UTRUE(idxq);
        % error
        e_h = uh - u;
        % store values
        xq_vec = [xq_vec; idxq];
        eh_vec = [eh_vec; e_h];
        uh_vec = [uh_vec; uh];
        u_vec  = [u_vec; u];
    end
end
gpar.xq = xq_vec;
gpar.eh_vec = eh_vec;
gpar.u_vec = u_vec;
gpar.uh_vec = uh_vec;
%%
if debug == 1
    % for testing uh function at the nodes
    uh = zeros(length(x),1);
    uexact = zeros(length(x),1);
    for j  = 1:N+1
        uh(j)  = UhFunc(x(j),j,x,h,c);
        uexact(j) = UTRUE(x(j));
    end
    uh(N+2)  = UhFunc(x(N+2),N+1,x,h,c);
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