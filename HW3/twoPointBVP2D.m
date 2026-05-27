function [errL2,errH1,MatAssemTime,solveMat] = twoPointBVP2D(par)

% define functions
if par.testfunc == 1
    fun.UTRUE   = @(x,y)  1 + 0*x;
    fun.UTRUEX  = @ (x,y) 0;
    fun.UTRUEXX = @ (x,y) 0;
    fun.UTRUEY  = @ (x,y) 0;
    fun.UTRUEYY = @ (x,y) 0;
elseif par.testfunc == 2
    fun.UTRUE   = @(x,y)  x;
    fun.UTRUEX  = @ (x,y) 1 + 0*x;
    fun.UTRUEXX = @ (x,y) 0;
    fun.UTRUEY  = @ (x,y) 0;
    fun.UTRUEYY = @ (x,y) 0;
elseif par.testfunc == 3
    fun.UTRUE   = @(x,y)  y;
    fun.UTRUEX  = @ (x,y) 0;
    fun.UTRUEXX = @ (x,y) 0;
    fun.UTRUEY  = @ (x,y) 1 + 0*y;
    fun.UTRUEYY = @ (x,y) 0;
else
    fun.UTRUE   = @ (x,y)  y.^3 + sin(5*(x+y)) + 2*exp(x);
    fun.UTRUEX  = @ (x,y) 5*cos(5*(x+y)) + 2*exp(x);
    fun.UTRUEXX = @ (x,y) -25*sin(5*(x+y)) + 2*exp(x);
    fun.UTRUEY  = @ (x,y) 3*y.^2 + 5*cos(5*(x+y));
    fun.UTRUEYY = @ (x,y) 6*y -25*sin(5*(x+y));
end
fun.FORCE = @(x,y) -par.p*(fun.UTRUEXX(x,y) + fun.UTRUEYY(x,y))...
    + par.q * (fun.UTRUE(x,y));

% mesh generation
[P, T] = Meshgen(par);
x = P(:,1);
y = P(:,2);
%% boundary and interior points
tol = 1e-12;
boundary_nodes = find(abs(P(:,1)) < tol | abs(P(:,1)-1) < tol | ...
    abs(P(:,2)) < tol | abs(P(:,2)-1) < tol );
interior_nodes = setdiff((1:size(P,1))',boundary_nodes);

%% form matrix and load vector
[A,F,Area,MatAssemTime] = formMatrixGlb(P,T,par,fun);

% Reduce matrix and load vector
AI = A(interior_nodes,interior_nodes);
B = A(interior_nodes,boundary_nodes);
Cb = fun.UTRUE(x(boundary_nodes),y(boundary_nodes));
FI = F(interior_nodes) - B * Cb;
tic
CI = AI\FI;
solveMat = toc;

%% solve for c
C = zeros(size(P,1),1);
C(boundary_nodes) = Cb;
C(interior_nodes) = CI;

%% compute error
[errL2,errH1] = getErrors(P,T,C,fun,Area);
%quadrature points
return
end

