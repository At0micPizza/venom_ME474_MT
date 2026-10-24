%% Q21: UD vs QUICK, Crank-Nicolson, dt = 0.0025, Ny = 20
clear; close all; clc

p = default_params();
p.theta = 0.5;  p.dt = 0.0025;  p.Ny = 20;  p.verbose = false;

Nxs     = [100 200 400 700];
schemes = {'UD', 'QUICK'};
outs    = cell(2, numel(Nxs));
pk      = zeros(2, numel(Nxs));
for s = 1:2
    for k = 1:numel(Nxs)
        p.scheme = schemes{s};  p.Nx = Nxs(k);
        p.name = sprintf('%s (%d,20)', schemes{s}, Nxs(k));
        outs{s,k} = venom_solver(p);
        pk(s,k) = max(outs{s,k}.phiOut);   % peak outlet concentration
    end
    fig = plot_results('timeseries', outs(s,:));
    plot_results('save', fig, ['q21_' lower(schemes{s})]);
end

% relative variation of the peak w.r.t. the finest mesh
relVar = abs(pk - pk(:,end))./pk(:,end);
fig = figure;
semilogy(Nxs(1:end-1), 100*relVar(1,1:end-1), 'o-', Nxs(1:end-1), 100*relVar(2,1:end-1), 's-');
grid on; xlabel('N_x'); ylabel('relative variation of max \phi_{out} [%]');
legend(schemes, 'Location', 'best');
plot_results('save', fig, 'q21_relative_variation');
