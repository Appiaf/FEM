function A = formMatrix(x,h,N,coeffOpt)
% Define Matrix A
A = zeros(N,N);
if coeffOpt == 1
    % First row
    A(1,1) = 3*(1/h(1) + 1/h(2)) + 2*(1/3*(h(1)+h(2)));
    A(1,2) = -3/h(2) + 2*(1/6*h(2));

    % Interior rows
    for i = 2:N-1
        A(i,i-1) = -3/h(i) + 1/3*h(i);
        A(i,i)   = 3*(1/h(i) + 1/h(i+1)) +2/3*(h(i)+h(i+1));
        A(i,i+1) = -3/h(i+1) + 1/3*(h(i+1));
    end

    % Last row
    A(N,N-1) = -3/h(N) + 1/3*(h(N));
    A(N,N)   = 3*(1/h(N) + 1/h(N+1)) + 2/3*(h(N)+h(N+1));
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