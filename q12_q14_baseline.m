%% Q12-Q14  Baseline: explicit Euler, (Nx,Ny) = (100,10), dt = 0.005 s, Pe = 10
clear; close all; clc

p = default_params();
p.name          = 'explicit, dt = 0.005 s';
p.snapshotTimes = [0.05 0.1 1 3 6 9 12];   % pick times that show each BC (Q13)
p.frameEvery    = 10;                      % one frame every 0.05 s
% note: after the bite phi drops to ~0.05, so a fixed [0 1] colour range hides
% the pulse. Keep ONE fixed range for all frames, but choose it sensibly.

out = venom_solver(p);

%% Q12: 2D field at a few times
fig = plot_results('snapshots', out, [0 1]);
plot_results('save', fig, 'q12_snapshots', 'image');

%% Q13: a-posteriori check of the boundary conditions
% TODO: inlet  -> phi(x=0,y) during and after the bite
%       outlet -> d(phi)/dx at x = L  (e.g. phi in the last 2-3 columns)
%       walls  -> phi(y) profiles near y = 0 and y = H, d(phi)/dy at the wall

%% Q14: outlet concentration and total venom
fig = plot_results('timeseries', {out}, 'Q14: baseline');
plot_results('save', fig, 'q14_timeseries');

% TODO: criterion for "venom has left the domain" and the resulting time,
%       e.g. first t after the peak where phiTot < 1% of max(phiTot)

%% video (check that the animation works; the real deliverable is Q24)
% plot_results('video', out, fullfile('output', 'q12_baseline.mp4'), 25, [0 0.06]);
