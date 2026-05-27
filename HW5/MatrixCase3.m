function [A0, A1, A2, M,K, F] = MatrixCase3(N, x)

A0 = sparse(N,N);
A1 = sparse(N,N);
A2 = sparse(N,N);
M = sparse(N,N);
K = sparse(N,N);
F  = zeros(N,1);

for k = 1:N+1
    h = x(k+1) - x(k);
    % Local stiffness
    Ke = (1/h)*[1 -1; -1 1];
    % Local mass
    Me  = h/6*[2 1;1 2];
    % Local advection
    Ade = 1/2*[-1 1;-1 1];
    
    % Load vector (f=1)
    Fe = (h/2)*[1;1];

    % Assembly
    nodes = [k-1, k];
    for ai = 1:2         %interior indexing
        i = nodes(ai);
        for bj = 1:2
            j = nodes(bj);
            if i >=1 && j >=1 && i <=N && j <=N
                A0(i,j) = A0(i,j) + Ke(ai,bj);
                A1(i,j) = A1(i,j) + Ade(ai,bj);
                A2(i,j) = A2(i,j) + Me(ai,bj);
                M(i,j)  = M(i,j)  + Me(ai,bj);
                K(i,j)  = K(i,j)  + Ke(ai,bj);
            end
        end
        if i >=1 && i <=N
            F(i) = F(i) + Fe(ai);
        end
    end
end

end