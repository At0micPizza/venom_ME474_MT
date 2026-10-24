function out = venom_solver(p)
% Unsteady 2D convection-diffusion of venom (FVM + theta method)
% usage: out = venom_solver(default_params())

if isempty(p.Gam)
    p.Gam = p.rho*p.umean*p.H/p.Pe;   % Gamma from the Peclet number (Q11)
end
m      = build_mesh(p);
nSteps = round(p.tEnd/p.dt);
t      = (0:nSteps)*p.dt;
phi    = p.phi0*ones(m.N, 1);         % initial condition

% time series
phiOut = nan(1, nSteps+1);
phiTot = nan(1, nSteps+1);
[phiOut(1), phiTot(1)] = compute_diagnostics(phi, m, p);

% stored 2D fields (snapshots for the plots, frames for the video)
snapK   = unique(min(max(round(p.snapshotTimes/p.dt), 0), nSteps));
snapPhi = nan(m.Ny, m.Nx, numel(snapK));
frameK  = [];
if p.frameEvery > 0
    frameK = 0:p.frameEvery:nSteps;
end
frames = nan(m.Ny, m.Nx, numel(frameK), 'single');
[snapPhi, frames] = store_fields(0, phi, m, snapK, snapPhi, frameK, frames);

% theta method: (M + theta*K) phi^{k+1} = (M - (1-theta)*K) phi^k
pulsating = p.f > 0;
[K, M] = assemble_operator(m, p, t(1));
if ~pulsating   % steady flow -> same matrices for every step
    A = apply_bc_matrix(M + p.theta*K, m, p);
    A = decomposition(A);   % factorize once, A\b then reuses the LU factors
    B = M - (1 - p.theta)*K;
end

if p.verbose
    fprintf('%s: theta = %g, dt = %g s, (Nx,Ny) = (%d,%d), Pe = %g, %d steps\n', ...
            p.name, p.theta, p.dt, m.Nx, m.Ny, p.Pe, nSteps);
end
diverged = false;
tic
for k = 1:nSteps
    if pulsating   % velocity changes in time -> rebuild K at t^{k+1}, keep the old one for t^k
        Kold   = K;
        [K, M] = assemble_operator(m, p, t(k+1));
        A = apply_bc_matrix(M + p.theta*K, m, p);
        B = M - (1 - p.theta)*Kold;
    end

    b   = apply_bc_rhs(B*phi, m, p, t(k+1));
    phi = A\b;

    [phiOut(k+1), phiTot(k+1)] = compute_diagnostics(phi, m, p);
    [snapPhi, frames] = store_fields(k, phi, m, snapK, snapPhi, frameK, frames);

    if any(~isfinite(phi)) || max(abs(phi)) > p.divergeLimit   % stop blown-up runs
        diverged = true;
        if p.verbose
            fprintf('  diverged at t = %.4g s\n', t(k+1));
        end
        break
    end
end
cpuTime = toc;
if p.verbose
    fprintf('  done in %.1f s\n', cpuTime);
end

out.p        = p;
out.mesh     = m;
out.t        = t;
out.phiOut   = phiOut;
out.phiTot   = phiTot;
out.phi      = phi;
out.snapT    = snapK*p.dt;
out.snapPhi  = snapPhi;
out.frameT   = frameK*p.dt;
out.frames   = frames;
out.diverged = diverged;
out.cpuTime  = cpuTime;

% Q14: venom has left once phiTot drops below 1% of its peak
[pk, kPk]  = max(phiTot);
kLeave     = find(phiTot(kPk:end) < 0.01*pk, 1) + kPk - 1;
out.tLeave = NaN;
if ~isempty(kLeave)
    out.tLeave = t(kLeave);
end
end


function m = build_mesh(p)
% nodes on all four boundaries, i = row (y), j = column (x), n = j + (i-1)*Nx  (Q6, Q7)
m.Nx = p.Nx;
m.Ny = p.Ny;
m.N  = p.Nx*p.Ny;
m.dx = p.L/(p.Nx-1);
m.dy = p.H/(p.Ny-1);
m.x  = (0:p.Nx-1)*m.dx;
m.y  = (0:p.Ny-1)'*m.dy;

m.dof = @(i,j) j + (i-1)*p.Nx;

% row and column index of every DOF
[J, I] = meshgrid(1:p.Nx, 1:p.Ny);
m.I = reshape(I', [], 1);
m.J = reshape(J', [], 1);

% DOFs touching each boundary (corners are in two lists)
m.bottom = m.dof(1, 1:p.Nx)';
m.top    = m.dof(p.Ny, 1:p.Nx)';
m.inlet  = m.dof(1:p.Ny, 1)';
m.outlet = m.dof(1:p.Ny, p.Nx)';
m.inner  = find(m.I > 1 & m.I < p.Ny & m.J > 1 & m.J < p.Nx);

% CV volumes per unit depth: half on the boundaries, quarter in the corners
wx = ones(1, p.Nx);  wx([1 end]) = 0.5;
wy = ones(p.Ny, 1);  wy([1 end]) = 0.5;
m.vol = reshape(((wy*wx)*m.dx*m.dy)', [], 1);

m.toField = @(phi) reshape(phi, p.Nx, p.Ny)';   % DOF vector -> Ny x Nx field
end


function ux = velocity_profile(y, t, p)
ux = 6*p.umean*(y/p.H).*(1 - y/p.H);   % Poiseuille profile (Q1)
if p.f > 0
    ux = max(p.alpha*ux*cos(2*pi*p.f*t), 0);   % pulsating flow (Q23)
end
end


function phi_in = inlet_concentration(t, p)
phi_in = p.phiBite*(t <= p.tBite + 1e-12);   % venom only during the bite, small tol for round-off
end


function [K, M] = assemble_operator(m, p, t)
% K*phi = a_P*phi_P - sum(a_nb*phi_nb) for the inner CVs, M = a_P0 on the diagonal (Q7)
ux  = velocity_profile(m.y, t, p);
Dx  = p.Gam*m.dy/m.dx;       % diffusive conductances (CD)
Dy  = p.Gam*m.dx/m.dy;
aP0 = p.rho*m.dx*m.dy/p.dt;  % unsteady term

n  = m.inner;
F  = p.rho*ux(m.I(n))*m.dy;  % mass flux through e and w faces, F_e = F_w
aE = Dx*ones(size(n));       % UD: u_x >= 0, so nothing is convected from E
aW = Dx + F;                 % UD: phi_w = phi_W
aWW = zeros(size(n));
aN = Dy*ones(size(n));
aS = Dy*ones(size(n));

if strcmp(p.scheme, 'QUICK')   % Q21, UD stays in column j = 2 (no WW node)
    q = m.J(n) >= 3;
    aE(q)  = Dx - 3/8*F(q);
    aW(q)  = Dx + 7/8*F(q);
    aWW(q) = -1/8*F(q);
end
aP = aE + aW + aWW + aN + aS;   % + (F_e - F_w) = 0

rows = [n; n; n; n; n; n];
cols = [n; n+1; n-1; n-2; n+m.Nx; n-m.Nx];
vals = [aP; -aE; -aW; -aWW; -aN; -aS];
keep = vals ~= 0;              % drops the WW entries of UD cells
K = sparse(rows(keep), cols(keep), vals(keep), m.N, m.N);
M = sparse(n, n, aP0, m.N, m.N);
end


function A = apply_bc_matrix(A, m, p)
% boundary rows, the inner rows are already in A
n = m.inlet;                       % inlet incl. corners: phi = phi_inlet (Q9)
A = A + sparse(n, n, 1, m.N, m.N);

n = m.outlet;                      % outlet incl. corners: phi_n - phi_{n-1} = 0 (Q8)
A = A + sparse([n; n], [n; n-1], [ones(size(n)); -ones(size(n))], m.N, m.N);

n = m.bottom(2:end-1);             % bottom wall: phi_n - phi_{n+Nx} = 0 (Q10)
A = A + sparse([n; n], [n; n+m.Nx], [ones(size(n)); -ones(size(n))], m.N, m.N);

n = m.top(2:end-1);                % top wall: phi_n - phi_{n-Nx} = 0 (Q10)
A = A + sparse([n; n], [n; n-m.Nx], [ones(size(n)); -ones(size(n))], m.N, m.N);
end


function b = apply_bc_rhs(b, m, p, t)
b(m.inlet)            = inlet_concentration(t, p);   % inlet value at t^{k+1} (Q9)
b(m.outlet)           = 0;                           % Neumann rows (Q8, Q10)
b(m.bottom(2:end-1))  = 0;
b(m.top(2:end-1))     = 0;
end


function [phiOut, phiTot] = compute_diagnostics(phi, m, p)
phiOut = trapz(m.y, phi(m.outlet))/p.H;   % y-averaged outlet concentration (Q14)
phiTot = sum(phi.*m.vol);                 % total venom per unit depth (Q14)
end


function [snapPhi, frames] = store_fields(k, phi, m, snapK, snapPhi, frameK, frames)
s = find(snapK == k, 1);
if ~isempty(s)
    snapPhi(:,:,s) = m.toField(phi);
end
f = find(frameK == k, 1);
if ~isempty(f)
    frames(:,:,f) = m.toField(phi);
end
end
