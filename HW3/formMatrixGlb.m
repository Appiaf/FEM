function [K, F,Area,MatAssemTime] = formMatrixGlb(P,T,par,fun)
q = par.q;
p = par.p;
N = size(P,1);
K = zeros(N,N);
F = zeros(N,1);
Area = zeros(size(T,1),1);
tic
for e = 1:size(T,1)

    nodes = T(e,:);
    Pe = P(nodes,:);

    [Ke,fe,areak] = formMatrixEle(Pe,p,q,fun);
    Area(e) = areak;
    % assemble
    for i = 1:3
        for j = 1:3
            K(nodes(i), nodes(j)) = ...
                K(nodes(i), nodes(j)) + Ke(i,j);
        end
    end
% contribution equally to nodes (P1 basis)
for i = 1:3
    F(nodes(i)) = F(nodes(i)) + fe(i);
end
end
MatAssemTime = toc;
return
end