function U_hx = UhFuncx(j,h,cj,N)
if j == 1
    U_hx = cj(1)/h(1);
elseif j == N+1
    U_hx = -cj(N)/h(N+1);
else
    U_hx = (-cj(j-1)+ cj(j))/h(j);
end
return
end