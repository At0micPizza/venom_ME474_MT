%% Q23-Q24: pulsating flow, absorbing walls, CN, dt = 0.0025 s, (700,20)
clear; close all; clc
p = default_params();
p.theta = 0.5;  p.dt = 0.0025;
p.Nx = 700;  p.Ny = 20;
p.snapshotTimes = [1 2 5];
p.C = 1e-4;

% run 1: steadz flow
p.f=0;
p.name = 'steady';
outS = venom_solver(p);

% run 2: pulsating flow

p.f = 1;
p.alpha = pi;          % mean of max(alpha*cos, 0) is alpha/pi -> same mean flow as steady
p.name = 'pulsating, f = 1Hz';
outP = venom_solver(p);

t = linspace(0,2,2001);
uc = 1.5*p.umean;
w = 2*pi*p.f;
vc = max(p.alpha*uc*cos(w*t),0);
fig = figure('Position', [100 100 600 300]); hold on; grid on; box on
plot(t, vc);
xlabel('t [s]');  ylabel('v_x(H/2, t) [m/s]');  title('Centreline velocity, f = 1 Hz');
plot_results('save', fig, 'q23_centreline');



%% step-14 plots: both runs
fig = plot_results('timeseries', {outS, outP});
lgd = findobj(fig, 'Type', 'Legend');  lgd.Location = 'southwest';
plot_results('save', fig, 'q23_timeseries');
fprintf('t_leave: steady %.2f s, pulsating %.2f s\n', outS.tLeave, outP.tLeave);



%% Q24
%24a - running heart rate f=3Hy

p.f = 3;
p.alpha = pi;
p.name = 'pulsating, f = 3Hz';
outP3 = venom_solver(p);


fig = plot_results('timeseries', {outS, outP, outP3});      % steady, 1 Hz, 3 Hz
lgd = findobj(fig, 'Type', 'Legend');  lgd.Location = 'northeast';
plot_results('save', fig, 'q24_timeseries_f');
fprintf('t_leave: steady %.2f s, 1 Hz %.2f s, 3 Hz %.2f s\n', outS.tLeave, outP.tLeave, outP3.tLeave);

% % Q24b - Pe = 100
% p.Pe = 100;
% p.frameEvery = 20;
% p.name = 'f=3Hz, Pe =100';
% outPe100 = venom_solver(p);
% 
% fig = plot_results('timeseries', {outP3, outPe100});        % Pe = 10 vs Pe = 100
% lgd = findobj(fig, 'Type', 'Legend');  lgd.Location = 'southwest';
% plot_results('save', fig, 'q24_timeseries_Pe');
% fprintf('t_leave: Pe = 10 %.2f s, Pe = 100 %.2f s\n', outP3.tLeave, outPe100.tLeave);
% 
% % Q24c video
% plot_results('video', outPe100, fullfile('output', 'q24_Pe100.mp4'), 25, [0 0.1]);