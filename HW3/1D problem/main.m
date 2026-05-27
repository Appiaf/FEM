paramSets = [
    1, 5;
    3, 10;
    5, 9
    ];
par.a = 0;  par.b = 1;
par.testOption = 1;
par.debug = 0;
%Nvalue = [9,19,39,79,159];
Nvalue = [10,20,40,80,160];
results = [];

for p = 1:size(paramSets,1)

    par.k1 = paramSets(p,1);
    par.k2 = paramSets(p,2);

    l2norm   = zeros(length(Nvalue),1);
    linfnorm = zeros(length(Nvalue),1);
    h = zeros(length(Nvalue),1);


    for j = 1:length(Nvalue)
        par.N = Nvalue(j);

        gpar = twoPointBVP(par);
        
        h(j) = 1/(Nvalue(j)+1);
        l2norm(j)   = gpar.L2Norm;
        linfnorm(j) = gpar.LinfNorm;
        %enorm(j)    = gpar.ENorm;

        if Nvalue(j) == 80
            figure(1);
            plot(gpar.xq, gpar.u_vec,'r-*','LineWidth', 1)
            hold on
            plot(gpar.xq, gpar.uh_vec,'b-o','LineWidth', 1)
            hold off
            grid on;
            xlim([0 1]);
            legend('Exact solution', 'FEM solution', 'Location','best')
            xlabel('$x$', 'Interpreter', 'latex');
            title(sprintf('Exact solution vs Approximate solution:$N+1=%d$, $k1=%g$, $k2=%g$', Nvalue(j)+1,par.k1,par.k2), ...
                'Interpreter','latex')
            filename = sprintf('img/bsolplot_N%d_k1%g_k2%g.eps',Nvalue(j)+1,par.k1,par.k2);
            print('-depsc2', filename);
        end
    end

    % convergence rates
    l2conv = NaN(length(Nvalue),1);
    linfconv = NaN(length(Nvalue),1);
    Econv = NaN(length(Nvalue),1);

    for j = 2:length(Nvalue)
        l2conv(j)   = log(l2norm(j-1)/l2norm(j)) ...
            / log((Nvalue(j)+1)/(Nvalue(j-1)+1));
        linfconv(j) = log(linfnorm(j-1)/linfnorm(j))...
            / log((Nvalue(j)+1)/(Nvalue(j-1)+1));
        % l2conv(j)   = log(l2norm(j-1)/l2norm(j))/log(2);
        % linfconv(j) = log(linfnorm(j-1)/linfnorm(j))/log(2);
        % Econv(j)    = log(enorm(j-1)/enorm(j))/log(2);
    end

    % store EVERYTHING (errors + rates)
    for j = 1:length(Nvalue)
        results = [results;
            par.k1, par.k2, Nvalue(j), ...
            l2norm(j), linfnorm(j), ...
            l2conv(j), linfconv(j)];
    end
end

T = array2table(results, ...
    'VariableNames', {'k1','k2','N', ...
    'L2Error','LinfError', ...
    'L2Conv','LinfConv'});

disp(T)