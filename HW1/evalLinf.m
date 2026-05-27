function Linfnorm = evalLinf(xval,h,cj,N)
%define functions
M = 3;
u   = @(x) x.*(x-1).*(sin(5*x)+3*exp(x));
eNorm = 0;
for j = 1: N+1
    for m = 1:M
        idx = xval(j) + h(j)*m/M;
        e_h = (UhFunc(idx,j,xval,h,cj,N) - u(idx));
        eNorm = max(eNorm ,abs(e_h));
    end
end
Linfnorm = eNorm;
return
end