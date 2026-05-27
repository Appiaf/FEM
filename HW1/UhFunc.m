function U_h = UhFunc(a,j,x,h,cj,N)
if j == 1
    U_h = 1/h(1)*(cj(1)*(a-x(1)));
elseif j == N+1
    U_h = (cj(N)*(x(N+2)-a))/h(N+1);
else
    U_h = (cj(j-1)*(x(j+1)-a)+cj(j)*(a-x(j)))/h(j);
end
return
end