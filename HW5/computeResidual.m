function r = computeResidual(mu, u_rom, A0, A1, A2, K, F,caseType)
if caseType == 1
    mu1 = mu(1); mu2 = mu(2);
    A = 2*A0 + mu1*A1 + mu2*A2;
elseif caseType == 2
    mu1 = mu(1); mu2 = mu(2);
    A = mu1*A0 + mu2* A1 + A2;
else
    A = 0.1*A0 + mu*A1 + A2;
end
res = F - A*u_rom;
r = sqrt(res'*(K\res));
end