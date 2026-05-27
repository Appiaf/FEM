function [A0_l, A1_r, A2, M,K, F] = MatrixCase2(N, x)

A0_l = sparse(N,N);
A1_r = sparse(N,N);
A2 = sparse(N,N);
M = sparse(N,N);
K = sparse(N,N);
F  = zeros(N,1);

for k = 1:N+1
    h = x(k+1) - x(k);
    xl = x(k); xr = x(k+1);
    % Local stiffness
    Ke = (1/h)*[1 -1; -1 1];
    % Local mass
    Mele  = h/6*[2 1;1 2];
    
    A0_le = zeros(2,2);
    A1_re = zeros(2,2);
    if xr <= 0.5
        % Local sin-weighted stiffness
        I = cos(x(k)) - cos(x(k+1));
        Ke1 = (1/h^2)*I * [1 -1; -1 1];
        A0_le = 2*Ke + Ke1;
    elseif xl >= 0.5
        I = 1/3*(x(k+1)^3 - x(k)^3);
        Ke1 = (1/h^2)*I * [1 -1; -1 1];
        A1_re = Ke1;
    else
        error('Choose mesh so x=0.5 is a node.');
    end

    % --- Reaction matrix using 2-point Gaussian quadrature ---
    Me = zeros(2,2);
    % Quadrature points (reference)
    xi = [-1/sqrt(3), 1/sqrt(3)];
    for q = 1:2
        % Map to physical element
        xq = (x(k) + x(k+1))/2 + (h/2)*xi(q);
        w = exp(xq);
        % Shape functions
        phi = [(1 - xi(q))/2; (1 + xi(q))/2];
        % Add contribution
        Me = Me + w * (phi * phi') * (h/2);
    end

    % Load vector (f=1)
    Fe = (h/2)*[1;1];

    % Assembly
    nodes = [k-1, k];
    for ai = 1:2         %interior indexing
        i = nodes(ai);
        for bj = 1:2
            j = nodes(bj);
            if i >=1 && j >=1 && i <=N && j <=N
                A0_l(i,j) = A0_l(i,j) + A0_le(ai,bj);
                A1_r(i,j) = A1_r(i,j) + A1_re(ai,bj);
                A2(i,j) = A2(i,j) + Me(ai,bj);
                M(i,j)  = M(i,j)  + Mele(ai,bj);
                K(i,j)  = K(i,j)  + Ke(ai,bj);
            end
        end
        if i >=1 && i <=N
            F(i) = F(i) + Fe(ai);
        end
    end
end

end