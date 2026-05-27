close all; clc;
caseType = 3;
basisType = 2;        % 1 Um  ; 2 QR
ErrEstimator = 1;     % 1 l1-norm   ; 2 residual
a = 0; b = 1;
xm = 0.5;
sigmaLB = 1.0;
N = 159;               % choose even so x=0.5 is a node
maxRB = 9;           % dim of Reduce Basis
nTest = 100;
tolGreedy = 1e-6;
x = linspace(a,b,N+2)';
h = (b-a)/(N+1);
%% Training set
n1 = 30; n2 = 30;
if caseType == 1
    sqL = 0; sqR = 1; sqB = 0; sqT = 2;
elseif caseType == 2
    sqL = 1; sqR = 4; sqB = 1; sqT = 4;
elseif caseType == 3
    sqL = 1;  sqR = 4;
end
mu1_vals = linspace(sqL,sqR,n1);
if caseType == 3
    Dtrain = mu1_vals(:);
    Dtest = sqL + sqR*rand(nTest,1);
else
    mu2_vals = linspace(sqB,sqT,n2);
    [mu1, mu2 ] = meshgrid(mu1_vals,mu2_vals);
    Dtrain = [mu1(:), mu2(:)];
    Dtest = [sqL + sqR*rand(nTest,1),sqB + sqT*rand(nTest,1)];
end
%% Assemble parameter-independent matrices
if caseType == 1
    [A0, A1, A2, M, K, F] = MatrixCase1(N, x);
elseif caseType == 2
    [A0, A1, A2, M, K, F] = MatrixCase2(N,x);
elseif caseType == 3
    [A0, A1, A2, M, K, F] = MatrixCase3(N,x);
end
%% INITIAL BASIS
if caseType == 3
    Initmu = 1/2*sqR;
else
    Initmu = [2,2];                 % row vector
end
Uh = solveFEM(Initmu, A0, A1, A2, F,caseType);
Um = Uh;                  % normalize
selected_mu = Initmu;             % store selected parameters
Delta_history = [];

%% GREEDY ALGORITHM
errH1 = zeros(maxRB+1,1);
errL2 = zeros(maxRB+1,1);
for m = 1:maxRB
    if basisType == 2
        [Qm, Rm] = qr(Um,0);
    else
        Qm = Um;
    end
    Delta = zeros(size(Dtrain,1),1);

    for i = 1:size(Dtrain,1)
        mu = Dtrain(i,:);

        % Skip already selected parameters
        if ismember(mu, selected_mu, 'rows')
            Delta(i) = -inf;
            continue;
        end

        % ROM solve
        [u_rom,c] = solveROM(mu, Qm, A0, A1, A2, F,caseType);

        if ErrEstimator == 1
            % L1 indicator
            if basisType == 1
                Delta(i) = norm(c,1);
            else
                Delta(i) = norm(Rm\c,1);
            end
        else
            % --- alternatively (better) ---
            r = computeResidual(mu, u_rom, A0, A1, A2, K, F,caseType);
            Delta(i) = r;
        end
    end

    % Find worst parameter
    [maxDelta, idx] = max(Delta);
    Delta_history = [Delta_history; maxDelta];

    fprintf('m = %d, max indicator = %.3e\n', m, maxDelta);

    if (maxDelta < tolGreedy) || (m == maxRB)
        % Compute error for current basis BEFORE breaking

        fprintf('\nComputing testing errors (final step)...\n');

        V = Qm;   % use current basis (already size m)

        maxErrH1 = 0;
        maxErrL2 = 0;

        for k = 1:nTest
            mutest = Dtest(k,:);
            UFOM = solveFEM(mutest, A0, A1, A2, F,caseType);
            UROM = solveROM(mutest, V, A0, A1, A2, F, caseType);

            e = UFOM - UROM;

            eH1 = sqrt(e' * K * e);
            eL2 = sqrt(e' * M  * e);

            maxErrH1 = max(maxErrH1, eH1);
            maxErrL2 = max(maxErrL2, eL2);
        end

        errH1(m) = maxErrH1;
        errL2(m) = maxErrL2;

        fprintf('m = %2d (final), H1 = %.4e, L2 = %.4e\n', ...
            m, errH1(m), errL2(m));

        break;
    end

    % New parameter
    mu_new = Dtrain(idx,:);

    % Store it
    selected_mu = [selected_mu; mu_new];

    % Full solve
    u_new = solveFEM(mu_new, A0, A1, A2, F,caseType);

    % Enrich basis
    Um = [Um, u_new];

    %% ----------------------------
    %  Testing errors
    %% ----------------------------
    fprintf('\nComputing testing errors...\n');

    V = Qm;
    %[Qm, ~] = qr(Um,0);

    maxErrH1 = 0;
    maxErrL2 = 0;

    for k = 1:nTest
        mutest = Dtest(k,:);
        UFOM = solveFEM(mutest, A0, A1, A2, F,caseType);

        UROM = solveROM(mutest, V, A0, A1, A2, F, caseType);

        e = UFOM - UROM;

        eH1 = sqrt(e' * K * e);
        eL2 = sqrt(e' * M  * e);

        maxErrH1 = max(maxErrH1, eH1);
        maxErrL2 = max(maxErrL2, eL2);
    end

    errH1(m) = maxErrH1;
    errL2(m) = maxErrL2;

    fprintf('m = %2d,   max H1-semi error = %.4e,   max L2 error = %.4e\n', ...
        m, errH1(m), errL2(m));
end
mStar = size(Um,2);
fprintf('\nOffline stage finished with m* = %d.\n', mStar);
%% Plot convergence
figure;
semilogy(1:length(Delta_history), Delta_history, '-o');
xlabel('m'); ylabel('max estimator/indicator');
title(sprintf('Greedy convergence: case = %d', caseType))
grid on;
filename = sprintf('img/GC_case%d_basis%d_err%d.eps',caseType,basisType,ErrEstimator);
print('-depsc2', filename);

figure;
semilogy(1:mStar, errH1(1:mStar), 'rs-','LineWidth',1.5); hold on;
semilogy(1:mStar, errL2(1:mStar), 'bo-','LineWidth',1.5);
grid on;
xlabel('m');
ylabel('max error');
title(sprintf('Testing errors: case = %d', caseType))
legend('semi-H^1 error','L^2 error','Location','best');
filename = sprintf('img/TE_case%d_basis%d_err%d.eps',caseType,basisType,ErrEstimator);
print('-depsc2', filename);

figure;
if caseType == 3
    plot(Dtrain, 0*Dtrain,'k.', 'MarkerSize', 8); hold on;
    plot(selected_mu, 0*selected_mu,'ro','LineWidth',1.5,'MarkerSize',7);
else
plot(Dtrain(:,1), Dtrain(:,2), 'k.', 'MarkerSize', 8); hold on;
plot(selected_mu(:,1), selected_mu(:,2), 'ro-','LineWidth',1.5,'MarkerSize',7);
end
grid on;
xlabel('\mu_1');
ylabel('\mu_2');
title(sprintf('Greedy selected parameter points: case =%d',caseType));
legend('Training set','Selected points','Location','best');
filename = sprintf('img/GSP_case%d_basis%d_err%d.eps',caseType,basisType,ErrEstimator);
print('-depsc2', filename);
%% ----------------------------
fprintf('\n================ SUMMARY ================\n');
fprintf('%d\n', caseType);
fprintf('Final RB dimension m*          = %d\n', mStar);
fprintf('Final greedy indicator         = %.4e\n', Delta_history(end));
fprintf('Final max test semi-H1 error   = %.4e\n', errH1(mStar));
fprintf('Final max test L2 error        = %.4e\n', errL2(mStar));
fprintf('=========================================\n');