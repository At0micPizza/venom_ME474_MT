%% Q16 + Q17: explicit vs implicit, implicit with different time steps
p = default_params();
p.name = 'explicit, dt = 0.005 s';
outEx = venom_solver(p);

p.theta = 1;
dts  = [0.01 0.005 0.0025 0.001 0.0005];
outs = cell(size(dts));
for k = 1:numel(dts)
    p.dt = dts(k);  p.name = sprintf('implicit, dt = %g s', dts(k));
    outs{k} = venom_solver(p);
end
fig = plot_results('timeseries', [{outEx}, outs]);
plot_results('save', fig, 'q16_q17_timeseries');

% Q18: Crank-Nicolson, dt = 0.0025
p = default_params();
p.theta = 0.5;  p.dt = 0.0025;  p.name = 'CN, dt = 0.0025 s';
outCN = venom_solver(p);

% reference with a very small dt, error = max difference of the outlet curve
p.dt = 0.0001;  p.name = 'CN, dt = 0.0001 s';
ref = venom_solver(p);
relErr = @(o) max(abs(o.phiOut - interp1(ref.t, ref.phiOut, o.t)))/max(ref.phiOut);
fprintf('Q18: CN dt = 0.0025 s: error = %.3f %%\n', 100*relErr(outCN));

% implicit Euler needs a much smaller dt for the same error
p.theta = 1;
dts = [0.0025 0.001 0.0005 0.00025 0.0001];
outsIE = cell(size(dts));
for k = 1:numel(dts)
    p.dt = dts(k);  p.name = sprintf('implicit, dt = %g s', dts(k));
    outsIE{k} = venom_solver(p);
    fprintf('     implicit dt = %g s: error = %.3f %%\n', dts(k), 100*relErr(outsIE{k}));
end
fig = plot_results('timeseries', [{outCN}, outsIE]);
plot_results('save', fig, 'q18_timeseries');
