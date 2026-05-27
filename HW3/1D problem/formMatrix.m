% function A = formMatrix(h,par,indk)
% % Define Matrix A
% k1 = par.k1;
% k2 = par.k2;
% N = par.N;
% A = zeros(N,N);
% % First row
% A(1,1) = k1*(1/h + 1/h);
% A(1,2) = -k1/h;
%
% % Interior rows
% for i = 2:N-1
%     if ( i < indk)
%         A(i,i-1) = -k1/h;
%         A(i,i)   = k1*(1/h + 1/h);
%         A(i,i+1) = -k1/h;
%     elseif i == indk
%         A(i,i-1) = -k1/h;
%         A(i,i)   = (k1/h + k2/h);
%         A(i,i+1) = -k2/h;
%     elseif (i > indk)
%         A(i,i-1) = -k2/h;
%         A(i,i)   = k2*(1/h + 1/h);
%         A(i,i+1) = -k2/h;
%     end
% end
%
% % Last row
% A(N,N-1) = -k2/h;
% A(N,N)   = k2*(1/h + 1/h);
% A = sparse(A);
% return
% end

function A = formMatrix(x,par)

k1 = par.k1;
k2 = par.k2;
N  = par.N;

A = zeros(N,N);
tol = 1e-12;

for e = 1:N+1   
    
    xL = x(e);
    xR = x(e+1);
    h  = xR - xL;
    
    % determine coefficient on element
    if xR <= 0.5 + tol
        % fully in left region
        coeff = k1 / h;
        
    elseif xL >= 0.5 - tol
        % fully in right region
        coeff = k2 / h;
        
    else
        % --- interface cuts this element ---
        
        h1 = 0.5 - xL;
        h2 = xR - 0.5;
        
        coeff = k1/h1 + k2/h2;  
        
    end
    
    % local stiffness matrix
    Ae = coeff * [1 -1; -1 1];
    
    % map to global indices (interior only)
    nodes = [e-1, e];   % interior indexing
    
    for a = 1:2
        for b = 1:2
            i = nodes(a);
            j = nodes(b);
            
            %interior points
            if i >= 1 && i <= N && j >= 1 && j <= N
                A(i,j) = A(i,j) + Ae(a,b);
            end
        end
    end
end

A = sparse(A);
end