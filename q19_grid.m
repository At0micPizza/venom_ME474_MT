%% Q19: grid refinement, Crank-Nicolson, dt = 0.0025
clear; close all; clc

p = default_params();
p.theta = 0.5;  p.dt = 0.0025;  p.verbose = false;

% y-direction
Nys = [10 20 30 40 50];
outsY = cell(size(Nys));
for k = 1:numel(Nys)
    p.Nx = 100;  p.Ny = Nys(k);  p.name = sprintf('(100,%d)', Nys(k));
    outsY{k} = venom_solver(p);
end

% x-direction
Nxs = 100:100:700;
outsX = cell(size(Nxs));
for k = 1:numel(Nxs)
    p.Nx = Nxs(k);  p.Ny = 20;  p.name = sprintf('(%d,20)', Nxs(k));
    outsX{k} = venom_solver(p);
end

% final check
p.Nx = 700;  p.Ny = 30;  p.name = '(700,30)';
out730 = venom_solver(p);
outsF = {outsX{end}, out730};

fprintf('%-10s %12s %12s %10s\n', 'mesh', 'max phiOut', 't of max', 'tLeave');
for o = [outsY, outsX, {out730}]
    [pk, k] = max(o{1}.phiOut);
    fprintf('%-10s %12.5f %12.3f %10.2f\n', o{1}.p.name, pk, o{1}.t(k), o{1}.tLeave);
end

fig = plot_results('timeseries', outsY);  plot_results('save', fig, 'q19_y_refinement');
fig = plot_results('timeseries', outsX);  plot_results('save', fig, 'q19_x_refinement');
fig = plot_results('timeseries', outsF);  plot_results('save', fig, 'q19_final_check');
