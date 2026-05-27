function A = formMatrix(x,h,par)
% Define Matrix A
coeffOpt = par.coeffOpt;
N = par.N;
A = zeros(N,N);
if coeffOpt == 1
    % First row
    A(1,1) = par.eps*(1/h + 1/h);
    A(1,2) = -par.eps/h + 1/2;

    % Interior rows
    for i = 2:N-1
        A(i,i-1) = -par.eps/h -1/2;
        A(i,i)   = par.eps*(1/h + 1/h);
        A(i,i+1) = -par.eps/h + 1/2;
    end

    % Last row
    A(N,N-1) = -par.eps/h -1/2;
    A(N,N)   = par.eps*(1/h + 1/h);
else
    % First row
    A(1,1) = (1+0.5*(x(1)+x(2)))/h(1) + (1+0.5*(x(2)+x(3)))/h(2);
    A(1,2) = -(1+0.5*(x(2)+x(3)))/h(2);

    % Interior rows
    for i = 2:N-1
        A(i,i-1) = -(1+0.5*(x(i)+x(i+1)))/h(i);
        A(i,i)   = (1+0.5*(x(i)+x(i+1)))/h(i) + (1+0.5*(x(i+1)+x(i+2)))/h(i+1);
        A(i,i+1) = -(1+0.5*(x(i+1)+x(i+2)))/h(i+1);
    end

    % Last row
    A(N,N-1) = -(1+0.5*(x(N+1)+x(N)))/h(N);
    A(N,N)   = (1+0.5*(x(N)+x(N+1)))/h(N) + (1+0.5*(x(N+1)+x(N+2)))/h(N+1);
end
A = sparse(A);
return
end