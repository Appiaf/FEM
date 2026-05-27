function lambda = nodalBasis(Pe,Ue,x)

lambda = zeros(3,1);
x1 = Pe(1,1); y1 = Pe(1,2);
x2 = Pe(2,1); y2 = Pe(2,2);
x3 = Pe(3,1); y3 = Pe(3,2);
xp = x(1); 
yp = x(2);

% area
area = 0.5 * det([x2-x1, x3-x1; y2-y1, y3-y1]);
% barycentric coordinates
lambda(1) = ((x2 - xp)*(y3 - yp) - (x3 - xp)*(y2 - yp)) / (2*area);
lambda(2) = ((xp - x1)*(y3 - y1) - (x3 - x1)*(yp - y1)) / (2*area);
lambda(3) = ((x2 - x1)*(yp - y1) - (xp - x1)*(y2 - y1)) / (2*area);

uh = Ue' * lambda;
return
end