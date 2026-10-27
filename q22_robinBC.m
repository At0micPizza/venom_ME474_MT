%% Q22: absorbing walls (Robin BC), CN, dt = 0.0025 s, (700,20)
clear; close all; clc
p = default_params();
p.theta = 0.5;  p.dt = 0.0025;
p.Nx = 700;  p.Ny = 20;
p.snapshotTimes = [1 2 5];

% run 1: closed walls
p.C = 0;     p.name = 'C = 0';
out0 = venom_solver(p);
% run 2: absorbing walls
p.C = 1e-4;  p.name = 'C = 1e-4';
out1 = venom_solver(p);
m = out1.mesh;

%% step-14 plots: both runs
fig = plot_results('timeseries', {out0, out1});
plot_results('save', fig, 'q22_timeseries');
fprintf('t_leave: C = 0: %.2f s, C = 1e-4: %.2f s\n', out0.tLeave, out1.tLeave);

%% profile across the channel at the pulse, t = 2 s
P0 = @(tt) out0.snapPhi(:,:,abs(out0.snapT - tt) < 1e-9);
P1 = @(tt) out1.snapPhi(:,:,abs(out1.snapT - tt) < 1e-9);
mid = round(m.Ny/2);
fig = figure('Position', [100 100 450 350]); hold on; grid on
F0 = P0(2);  F1 = P1(2);
[~, jPk] = max(F0(mid,:));                      % column at the pulse
plot(F0(:,jPk), m.y*1e3, 'o-', F1(:,jPk), m.y*1e3, 's-');
xlabel(sprintf('\\phi(x = %.2f m, y) [g/L]', m.x(jPk))); ylabel('y [mm]');
legend('C = 0', 'C = 1e-4');  title('walls, t = 2 s');
plot_results('save', fig, 'q22_profile');