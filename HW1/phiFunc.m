function phi = phiFunc(a,j,x,h)
xj   = x(j+1); xjm1 = x(j); xjp1 = x(j+2);
if (a >= xjm1) && (a < xj)
    phi = (a-xjm1)/h(j);
elseif (a >= xj) && (a <= xjp1)
    phi = (xjp1 - a)/h(j+1);
else
    phi = 0;
end
return
end