%% Q15 + Q20: explicit Euler with larger time steps / finer mesh
clear; close all; clc

% Q15: (100,10) with dt = 0.005 and dt = 0.01
p = default_params();
p.name = 'explicit, dt = 0.005 s';
out1 = venom_solver(p);
p.dt = 0.01;  p.name = 'explicit, dt = 0.01 s';
out2 = venom_solver(p);

fig = plot_results('timeseries', {out1, out2});
plot_results('save', fig, 'q15_timeseries');

% boundedness: a_P0 - a_P >= 0 for every inner CV
dtMax = @(p) 1/(2*p.Gam/(p.rho*(p.L/(p.Nx-1))^2) + 2*p.Gam/(p.rho*(p.H/(p.Ny-1))^2) ...
                + 1.5*p.umean/(p.L/(p.Nx-1)));
fprintf('Q15: dt_max = %.5f s on (%d,%d)\n', dtMax(out1.p), p.Nx, p.Ny);

% Q20: (700,20) with dt = 0.005
p = default_params();
p.Nx = 700;  p.Ny = 20;  p.name = 'explicit (700,20), dt = 0.005 s';
out3 = venom_solver(p);
fprintf('Q20: dt_max = %.5f s on (700,20), diverged = %d\n', dtMax(out3.p), out3.diverged);
