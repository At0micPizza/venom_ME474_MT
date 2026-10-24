%% Q11-Q14: baseline, explicit Euler, (Nx,Ny) = (100,10), dt = 0.005 s, Pe = 10
clear; close all; clc

p = default_params();

% Q11
Re  = p.H*p.umean/p.nu;
Gam = p.rho*p.umean*p.H./[10 100];
fprintf('Re = %.1f, Gamma = %.4f (Pe = 10) to %.4f (Pe = 100) kg/(m s)\n', Re, Gam);

% Q12
p.name = 'explicit, dt = 0.005 s';
p.snapshotTimes = [0.05 0.1 0.5 2 5 9 12];
out = venom_solver(p);
m = out.mesh;

fig = plot_results('snapshots', out, [0 0.1]);
plot_results('save', fig, 'q12_snapshots', 'image');

%% Q13: boundary conditions a posteriori
P = @(tt) out.snapPhi(:,:,abs(out.snapT - tt) < 1e-9);   % field at a stored time
mid = round(m.Ny/2);

fig = figure('Position', [100 100 900 320]);
subplot(1,3,1); hold on; grid on
Pa = P(0.05); Pb = P(0.5);
plot(Pa(:,1), m.y*1e3, 'o-', Pb(:,1), m.y*1e3, 's-');
xlabel('\phi(x = 0, y) [g/L]'); ylabel('y [mm]'); xlim([-0.1 1.1]);
legend('t = 0.05 s', 't = 0.5 s', 'Location', 'best'); title('inlet: \phi = \phi_{inlet}(t)');

subplot(1,3,2); hold on; grid on
Pc = P(9);
plot(m.x(end-10:end), Pc(mid, end-10:end), 'o-');
xlabel('x [m]'); ylabel('\phi(x, H/2) [g/L]'); title('outlet: \partial\phi/\partial x = 0, t = 9 s');

subplot(1,3,3); hold on; grid on
Pd = P(2);
[~, jPk] = max(Pd(mid,:));                      % column at the pulse
plot(Pd(:,jPk), m.y*1e3, 'o-');
xlabel(sprintf('\\phi(x = %.2f m, y) [g/L]', m.x(jPk))); ylabel('y [mm]');
title('walls: \partial\phi/\partial y = 0, t = 2 s');
plot_results('save', fig, 'q13_bc_check');

%% Q14
fig = plot_results('timeseries', {out});
plot_results('save', fig, 'q14_timeseries');
fprintf('venom left the domain (phiTot < 1%% of peak) at t = %.2f s\n', out.tLeave);
