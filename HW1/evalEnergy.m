function Energynorm = evalEnergy(xval,h,cj,N)
%define functions
M = 3;
ux   = @(x) x*(x-1)*(5*cos(5*x) + 3*exp(x)) ...
    + (sin(5*x)+3*exp(x))*(2*x-1);
eNorm = 0;
for j = 1: N+1
    for m = 1:M
        idx = xval(j) + m*h(j)/M;
        e_h = (UhFuncx(j,h,cj,N) - ux(idx))^2 * h(j)/M;
        eNorm = eNorm + e_h;
    end
end
Energynorm = sqrt(eNorm);
return
end