%ind : interior index
%P : vector of nodes
%T : matrix of elements
%area : vector of areas;
function [Ak,loadVecEle,K] = formMatrixEle(Pe,p,q,fun)
x1 = Pe(1,1); y1 = Pe(1,2);
x2 = Pe(2,1); y2 = Pe(2,2);
x3 = Pe(3,1); y3 = Pe(3,2);

K = 0.5 * abs((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));

% b(x cordinate for grad of lambda) and c (y cordinate) coefficients
b = [y2 - y3; y3 - y1; y1 - y2];
c = [x3 - x2; x1 - x3; x2 - x1];
%% Matrix components
Sk = zeros(3,3);
Mk = K/12 * [2 1 1;1 2 1;1 1 2];
for i = 1:3
    for j = 1:3
        Sk(i,j) = (b(i)*b(j) + c(i)*c(j)) / (4*K);
    end
end
Ak = q*Mk + p*Sk;
fe = [fun.FORCE(x1,y1);fun.FORCE(x2,y2);fun.FORCE(x3,y3)];
loadVecEle = Mk * fe;
return
end