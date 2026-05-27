function U_h = UhFunc(a,j,x,h,c)
U_h = (c(j)*(x(j+1)-a)+c(j+1)*(a-x(j)))/h(j);
return
end