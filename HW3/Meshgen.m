function[P,T] = Meshgen(par)
N = par.N;
xa = par.xa;
ya = par.ya;
meshType = par.meshType;
Nxtot = N + 2;
Nytot = N + 2;
Ntot  = Nytot*Nxtot;
Nele  = 2*(N+1)*(N+1);
P = zeros(Ntot,2);
T = zeros(Nele,3);
x = zeros(Nxtot,1);
y = zeros(Nytot,1);
x(1) = xa;
y(1) = ya;
ind = @(i,j) (j-1)*Nxtot + i;
h = 1/(N+1);
for j = 1 : N+1
    x(j+1) = h + x(j);
    y(j+1) = h + y(j);
end

for j = 1 : Nxtot
    for i = 1:Nytot
        ind_glb = ind(i,j);
        P(ind_glb,1) = x(i);
        P(ind_glb,2) = y(j);
    end
end
%% perturbed
if strcmp(meshType,'nonuniform')
    tol = 1e-12;
    for k = 1:size(P,1)

        x = P(k,1);
        y = P(k,2);

        % check if interior node
        if (x > tol && x < 1 - tol && y > tol && y < 1 - tol)

            Rx = rand - 0.5;
            Ry = rand - 0.5;

            P(k,1) = P(k,1) + (h/10)*Rx;
            P(k,2) = P(k,2) + (h/10)*Ry;

        end

    end
end
% elements
e = 1;
for j = 1:N + 1
    for i = 1:N + 1

        n1 = (j-1)*(Nxtot) + i;
        n2 = n1 + 1;
        n3 = n1 + Nxtot;
        n4 = n3 + 1;

        T(e,:) = [n1 n2 n3];
        T(e+1,:) = [n2 n4 n3];

        e = e+2;

    end
end
if par.N < 5
    triplot(T, P(:,1), P(:,2),'Linewidth',1.5)
    %axis equal
    hold on

    for i = 1:size(P,1)
        text(P(i,1), P(i,2), sprintf('%d', i), 'Color','r')
    end
hold off
title(sprintf('%s Mesh Grid; N = %d',meshType, N))
filename = sprintf('img/mesh_%s.eps',meshType);
print('-depsc2', filename);
end
return
end