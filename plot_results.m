function varargout = plot_results(mode, varargin)
% PLOT_RESULTS  All plotting and export in one place.
%
%   fig = plot_results('timeseries', {out1, out2, ...}, title)
%         Step-14 plot: y-averaged outlet concentration and total venom vs t.
%         Legend from out.p.name; diverged runs are marked and not used for the y-limits.
%
%   fig = plot_results('snapshots', out, clim)
%         All stored snapshots (p.snapshotTimes) stacked, same colour range.
%
%   h   = plot_results('field', ax, mesh, Phi, t, clim)
%         One filled contour plot of Phi (Ny x Nx). Aspect ratio NOT respected.
%
%         plot_results('video', out, filename, fps, clim)
%         mp4 of the stored frames (run with p.frameEvery > 0).
%
%         plot_results('save', fig, name)            vector PDF -> output/<name>.pdf
%         plot_results('save', fig, name, 'image')   raster PDF (use for contour fields)
%
%   clim is the fixed colour range [min max]; default [0 1].

switch lower(mode)
    case 'timeseries', varargout = {plot_timeseries(varargin{:})};
    case 'snapshots',  varargout = {plot_snapshots(varargin{:})};
    case 'field',      varargout = {plot_field(varargin{:})};
    case 'video',      make_video(varargin{:});
    case 'save',       save_figure(varargin{:});
    otherwise, error('plot_results:mode', 'Unknown mode "%s".', mode);
end
end


function fig = plot_timeseries(outs, ttl)
if ~iscell(outs), outs = {outs}; end

fig = figure('Position', [100 100 800 600]);
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
if nargin >= 2 && ~isempty(ttl), title(tl, ttl); end

ax1 = nexttile(tl); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
yl1 = [inf -inf];  yl2 = [inf -inf];
for r = 1:numel(outs)
    o = outs{r};
    lbl = o.p.name;
    if o.diverged
        lbl = [lbl ' (diverged)']; %#ok<AGROW>
    else
        yl1 = [min(yl1(1), min(o.phiOut)) max(yl1(2), max(o.phiOut))];
        yl2 = [min(yl2(1), min(o.phiTot)) max(yl2(2), max(o.phiTot))];
    end
    plot(ax1, o.t, o.phiOut, 'DisplayName', lbl);
    plot(ax2, o.t, o.phiTot, 'DisplayName', lbl);
end
if all(isfinite([yl1 yl2])) && yl1(2) > yl1(1) && yl2(2) > yl2(1)
    ylim(ax1, yl1 + [-0.05 0.05]*diff(yl1));
    ylim(ax2, yl2 + [-0.05 0.05]*diff(yl2));
end
ylabel(ax1, '$\bar\phi_{out}$ [g/L]', 'Interpreter', 'latex');
ylabel(ax2, '$\phi_{Total}$ [g/L$\cdot$m$^2$]', 'Interpreter', 'latex');
xlabel(ax2, 't [s]');
legend(ax1, 'Location', 'best', 'Interpreter', 'none');
linkaxes([ax1 ax2], 'x');
end


function fig = plot_snapshots(out, clim_)
if nargin < 2, clim_ = []; end
nS  = numel(out.snapT);
fig = figure('Position', [100 100 900 max(250, 160*nS)]);
tl  = tiledlayout(fig, nS, 1, 'TileSpacing', 'compact');
for s = 1:nS
    plot_field(nexttile(tl), out.mesh, out.snapPhi(:,:,s), out.snapT(s), clim_);
end
end


function h = plot_field(ax, m, Phi, tNow, clim_)
if nargin < 5 || isempty(clim_), clim_ = [0 1]; end

if max(Phi, [], 'all') > min(Phi, [], 'all')
    [~, h] = contourf(ax, m.x, m.y*1e3, Phi, 30, 'LineStyle', 'none');
else
    h = imagesc(ax, m.x, m.y*1e3, Phi);   % constant field: contourf cannot draw it
    axis(ax, 'xy');
end
clim(ax, clim_);
colormap(ax, turbo);
cb = colorbar(ax);
cb.Label.String = '\phi [g/L]';
xlabel(ax, 'x [m]');
ylabel(ax, 'y [mm]');
if nargin >= 4 && ~isempty(tNow)
    title(ax, sprintf('t = %.3f s', tNow));
end
end


function make_video(out, filename, fps, clim_)
if nargin < 3 || isempty(fps),   fps   = 25;    end
if nargin < 4 || isempty(clim_), clim_ = [0 1]; end
if isempty(out.frameT)
    error('plot_results:noFrames', 'No frames stored: run the solver with p.frameEvery > 0.');
end

vw = VideoWriter(filename, 'MPEG-4');
vw.FrameRate = fps;
open(vw);
cleanup = onCleanup(@() close(vw));

fig = figure('Position', [100 100 1000 300], 'Color', 'w', 'Visible', 'off');
ax  = axes(fig);
for f = 1:numel(out.frameT)
    if all(isnan(out.frames(:,:,f)), 'all'), break; end   % run diverged / stopped
    cla(ax);
    plot_field(ax, out.mesh, double(out.frames(:,:,f)), out.frameT(f), clim_);
    writeVideo(vw, getframe(fig));
end
close(fig);
fprintf('video written: %s\n', filename);
end


function save_figure(fig, name, content)
if nargin < 3, content = 'vector'; end
outDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~isfolder(outDir), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir, [name '.pdf']), 'ContentType', content, 'Resolution', 300);
end
