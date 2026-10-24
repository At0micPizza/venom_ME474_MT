%% Q19 + Q21: grid refinement with UD and QUICK, Crank-Nicolson, dt = 0.0025
clear; close all; clc

p = default_params();
p.theta = 0.5;  p.dt = 0.0025;  p.verbose = false;

% Q19: y-direction, UD
Nys = [10 20 30 40 50];
outsY = cell(size(Nys));
for k = 1:numel(Nys)
    p.Nx = 100;  p.Ny = Nys(k);  p.name = sprintf('(100,%d)', Nys(k));
    outsY{k} = venom_solver(p);
end

% Q19 + Q21: x-direction, UD and QUICK
Nxs     = 100:100:700;
schemes = {'UD', 'QUICK'};
outsX   = cell(2, numel(Nxs));
for s = 1:2
    for k = 1:numel(Nxs)
        if s == 2 && ~ismember(Nxs(k), [100 200 400 700])   % QUICK only on the Q21 meshes
            continue
        end
        p.scheme = schemes{s};  p.Nx = Nxs(k);  p.Ny = 20;
        p.name = sprintf('(%d,20)', Nxs(k));
        outsX{s,k} = venom_solver(p);
    end
end

% Q19: final check
p.scheme = 'UD';  p.Nx = 700;  p.Ny = 30;  p.name = '(700,30)';
out730 = venom_solver(p);

fprintf('%-10s %-6s %12s %12s %10s\n', 'mesh', 'scheme', 'max phiOut', 't of max', 'tLeave');
for o = [outsY, outsX(~cellfun(@isempty, outsX))', {out730}]
    [pk, k] = max(o{1}.phiOut);
    fprintf('%-10s %-6s %12.5f %12.3f %10.2f\n', o{1}.p.name, o{1}.p.scheme, pk, o{1}.t(k), o{1}.tLeave);
end

% relative variation of the peak w.r.t. (700,20), on the Q21 meshes
k21 = [1 2 4];
pk  = cellfun(@(o) max(o.phiOut), outsX(:, [k21 7]));
relVar = abs(pk(:,1:3) - pk(:,4))./pk(:,4);

% one figure: outlet concentration for each study + relative variation
fig = figure('Position', [100 100 900 600]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
studies = {outsY, outsX(1,:), outsX(2,~cellfun(@isempty, outsX(2,:)))};
titles  = {'(a) UD, refinement in y', '(b) UD, refinement in x', '(c) QUICK, refinement in x'};
for s = 1:3
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    for o = studies{s}
        plot(ax, o{1}.t, o{1}.phiOut, 'DisplayName', o{1}.p.name);
    end
    xlim(ax, [6 14]); ylim(ax, [0 0.09]);
    xlabel(ax, 't [s]'); ylabel(ax, '$\bar\phi_{out}$ [g/L]', 'Interpreter', 'latex');
    title(ax, titles{s}); legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
end
ax = nexttile(tl);
semilogy(ax, Nxs(k21), 100*relVar(1,:), 'o-', Nxs(k21), 100*relVar(2,:), 's-');
grid(ax, 'on'); xlabel(ax, 'N_x'); ylabel(ax, 'rel. variation of max \phi_{out} [%]');
title(ax, '(d) w.r.t. (700,20)'); legend(ax, schemes, 'Location', 'best');
plot_results('save', fig, 'q19_q21_grid');
