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

%% venom removed by the walls: difference of phiTot between the two runs
fig = figure('Position', [100 100 600 300]); hold on; grid on; box on
plot(out0.t, out0.phiTot - out1.phiTot);                  % C = 0 minus C = 1e-4
xlabel('t [s]');
ylabel('$\phi_{Total}(C=0) - \phi_{Total}(C=10^{-4})$', 'Interpreter', 'latex');
title('Venom removed by the walls');
plot_results('save', fig, 'q22_absorbed');

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

%% Q22: relative effect of absorption
relEffect = max(abs(out0.phiTot - out1.phiTot)) / max(out0.phiTot);
fprintf('relative change of phiTot due to absorption: %.2e\n', relEffect);

%% Q22: a posteriori check of the Robin BC at t = 2 s
Phi = P1(2);                                       % absorbing run, stored field at t = 2 s
Gam = out1.p.Gam;   C = out1.p.C;
lhs_bot = Gam * (Phi(2,:)     - Phi(1,:))   / m.dy;   % Gamma*(phi_N - phi_P)/dy, bottom wall
lhs_top = Gam * (Phi(end-1,:) - Phi(end,:)) / m.dy;   % Gamma*(phi_S - phi_P)/dy, top wall
fig = figure('Position', [100 100 700 350]); hold on; grid on; box on
plot(m.x, lhs_bot, '-',  m.x, C*Phi(1,:),   '--');
plot(m.x, lhs_top, '-',  m.x, C*Phi(end,:), ':');
xlim([0 0.5]);  xlabel('x [m]');  ylabel('\Gamma \partial\phi/\partial n  and  C\phi');
legend('bottom: \Gamma (\phi_N-\phi_P)/\Delta y', 'bottom: C\phi_P', ...
    'top: \Gamma (\phi_S-\phi_P)/\Delta y',    'top: C\phi_P', 'Location', 'best');
title('Robin BC a posteriori, t = 2 s');
plot_results('save', fig, 'q22_bc_check');