function [A0, A1, A2, M, K, F] = assemble_matrices(N, x)

A0 = sparse(N,N);
A1 = sparse(N,N);
A2 = sparse(N,N);
M = sparse(N,N);
F  = zeros(N,1);

for k = 1:N+1

    h = x(k+1) - x(k);

    % Local stiffness
    Ke = (1/h)*[1 -1; -1 1];

    % Local sin-weighted stiffness
    %I = integral(@(xx) sin(xx), x(k), x(k+1));
    I = cos(x(k)) - cos(x(k+1));
    Ke1 = (1/h^2)*I * [1 -1; -1 1];

    Mele  = h/6*[2 1;1 2];

    % Local mass (reaction)
    % --- Reaction matrix using 2-point Gaussian quadrature ---
    Me = zeros(2,2);
    % Quadrature points (reference)
    xi = [-1/sqrt(3), 1/sqrt(3)];

    for q = 1:2

        % Map to physical element
        xq = (x(k) + x(k+1))/2 + (h/2)*xi(q);

        % Only integrate on [0, 0.5]
        if xq <= 0.5 + 1e-12

            w = exp(xq);

            % Shape functions
            phi = [(1 - xi(q))/2; (1 + xi(q))/2];

            % Add contribution
            Me = Me + w * (phi * phi') * (h/2);
        end
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
                A0(i,j) = A0(i,j) + Ke(ai,bj);
                A1(i,j) = A1(i,j) + Ke1(ai,bj);
                A2(i,j) = A2(i,j) + Me(ai,bj);
                M(i,j)  = M(i,j)  + Mele(ai,bj);
            end
        end
        if i >=1 && i <=N
            F(i) = F(i) + Fe(ai);
        end
    end
end
K = A0;
end