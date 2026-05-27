function phi_dev = devphiFunc(a,j,x,h)
xj   = x(j+1); xjm1 = x(j); xjp1 = x(j+2);
if (a >= xjm1) && (a < xj)
    phi_dev = (1)/h(j);
elseif (a >= xj) && (a <= xjp1)
    phi_dev = (-1)/h(j+1);
else
    phi_dev = 0;
end
return
end