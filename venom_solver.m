function out = venom_solver(p)
% VENOM_SOLVER  Unsteady 2D convection-diffusion of venom in a plane channel.
%   out = venom_solver(p)      with p from default_params()
%
%   Theta scheme (report Q4-Q5) for every row that carries a finite-volume balance:
%       (M + theta*K^{k+1}) phi^{k+1} = (M - (1-theta)*K^k) phi^k
%   Rows of boundary CVs are then replaced by their boundary equation
%   (apply_bc_matrix for A, apply_bc_rhs for b) and A*phi = b is solved with backslash.
%
%   out fields:
%       t, phiOut, phiTot   time series (y-averaged outlet concentration, total venom)
%       phi                 final DOF vector
%       snapT, snapPhi      2D fields (Ny x Nx x nSnap) at p.snapshotTimes
%       frameT, frames      2D fields every p.frameEvery steps (for the video)
%       mesh, p, diverged, cpuTime
%
%   Local functions in this file:
%       build_mesh, velocity_profile, inlet_concentration     (done)
%       assemble_operator, apply_bc_matrix, apply_bc_rhs,     (TODO)
%       compute_diagnostics                                   (TODO)

if isempty(p.Gam)
    p.Gam = p.rho*p.umean*p.H/p.Pe;
end
m      = build_mesh(p);
nSteps = round(p.tEnd/p.dt);
t      = (0:nSteps)*p.dt;

phi = p.phi0*ones(m.N, 1);

%% storage
phiOut = nan(1, nSteps+1);
phiTot = nan(1, nSteps+1);
[phiOut(1), phiTot(1)] = compute_diagnostics(phi, m, p);

snapK   = unique(min(max(round(p.snapshotTimes/p.dt), 0), nSteps));
snapPhi = nan(m.Ny, m.Nx, numel(snapK));
if p.frameEvery > 0
    frameK = 0:p.frameEvery:nSteps;
else
    frameK = [];
end
frames = nan(m.Ny, m.Nx, numel(frameK), 'single');
[snapPhi, frames] = store_fields(0, phi, m, snapK, snapPhi, frameK, frames);

%% matrices
% With a steady flow K, M and the boundary rows do not change: build A once.
% With a pulsating flow K depends on t and is rebuilt every step.
pulsating = p.f > 0;
[K, M] = assemble_operator(m, p, t(1));
if ~pulsating
    A = apply_bc_matrix(M + p.theta*K, m, p);
    B = M - (1 - p.theta)*K;
end

%% time marching
if p.verbose
    fprintf('%s: theta = %g, dt = %g s, (Nx,Ny) = (%d,%d), Pe = %g, %d steps\n', ...
            p.name, p.theta, p.dt, m.Nx, m.Ny, p.Pe, nSteps);
end
diverged = false;
tic
for k = 1:nSteps
    if pulsating
        Kold   = K;                               % K at t^k
        [K, M] = assemble_operator(m, p, t(k+1)); % K at t^{k+1}
        A = apply_bc_matrix(M + p.theta*K, m, p);
        B = M - (1 - p.theta)*Kold;
    end

    b   = apply_bc_rhs(B*phi, m, p, t(k+1));
    phi = A\b;

    [phiOut(k+1), phiTot(k+1)] = compute_diagnostics(phi, m, p);
    [snapPhi, frames] = store_fields(k, phi, m, snapK, snapPhi, frameK, frames);

    if any(~isfinite(phi)) || max(abs(phi)) > p.divergeLimit
        diverged = true;
        if p.verbose
            fprintf('  diverged at t = %.4g s (max|phi| = %.3g)\n', t(k+1), max(abs(phi)));
        end
        break
    end
end
cpuTime = toc;
if p.verbose
    fprintf('  done in %.1f s\n', cpuTime);
end

%% output
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
end


%% =========================================================================
%  Mesh, velocity, inlet
%  =========================================================================

function m = build_mesh(p)
% Uniform structured mesh with nodes ON all four boundaries.
%   Indices (report Q6):  i = row    (y),  i = 1 bottom wall ... Ny top wall
%                         j = column (x),  j = 1 inlet       ... Nx outlet
%   DOF numbering (Q7), row by row:   n = j + (i-1)*Nx
%       neighbours of n:  E = n+1,  W = n-1,  N = n+Nx,  S = n-Nx

m.Nx = p.Nx;  m.Ny = p.Ny;  m.N = p.Nx*p.Ny;
m.dx = p.L/(p.Nx-1);
m.dy = p.H/(p.Ny-1);
m.x  = (0:p.Nx-1)*m.dx;            % 1 x Nx node positions
m.y  = (0:p.Ny-1)'*m.dy;           % Ny x 1 node positions

m.dof = @(i,j) j + (i-1)*p.Nx;     % (row, column) -> DOF index

% (i,j) of every DOF, as N x 1 vectors
[J, I] = meshgrid(1:p.Nx, 1:p.Ny);
m.I = reshape(I', [], 1);
m.J = reshape(J', [], 1);

% DOF lists of the CVs touching each boundary (corners appear in two lists)
m.bottom = m.dof(1,    1:p.Nx);
m.top    = m.dof(p.Ny, 1:p.Nx);
m.inlet  = m.dof(1:p.Ny, 1);
m.outlet = m.dof(1:p.Ny, p.Nx);
m.inner  = find(m.I > 1 & m.I < p.Ny & m.J > 1 & m.J < p.Nx);

% CV volume per unit depth: half CVs on the boundaries, quarter CVs in corners
wx = ones(1, p.Nx);  wx([1 end]) = 0.5;
wy = ones(p.Ny, 1);  wy([1 end]) = 0.5;
m.vol = reshape(((wy*wx)*m.dx*m.dy)', [], 1);   % sum(m.vol) == L*H

m.toField = @(phi) reshape(phi, p.Nx, p.Ny)';   % DOF vector -> Ny x Nx field Phi(i,j)
end


function ux = velocity_profile(y, t, p) %#ok<DEFNU>  (used by assemble_operator)
% Streamwise velocity u_x(y,t) [m/s]  (u_y = 0 everywhere).
%   Steady (p.f == 0): plane Poiseuille flow, Q1.
%   Pulsating (p.f > 0): v_x = max(alpha*u_x(y)*cos(2*pi*f*t), 0), Q23.

ux = 6*p.umean*(y/p.H).*(1 - y/p.H);

if p.f > 0
    if isempty(p.alpha)
        error('venom_solver:alpha', 'Set p.alpha for the pulsating flow (Q23).');
    end
    ux = max(p.alpha*ux*cos(2*pi*p.f*t), 0);
end
end


function phi_in = inlet_concentration(t, p) %#ok<DEFNU>  (used by apply_bc_rhs)
% Venom concentration imposed at the inlet [g/L]  (Q9, Q12).
phi_in = p.phiBite * (t <= p.tBite + 1e-12);   % tolerance: t = k*dt is not exact in floating point
end


%% =========================================================================
%  Discretization  (TODO)
%  =========================================================================

function [K, M] = assemble_operator(m, p, t) %#ok<INUSD,STOUT>
% Spatial operator K and "mass" matrix M of the theta scheme.
%
%   Each CV with a finite-volume balance satisfies (report Q5)
%       M*(phi^{k+1} - phi^k) + theta*K*phi^{k+1} + (1-theta)*K*phi^k = 0
%   so for row n:
%       K(n,n)  =  a_P           K(n,nb) = -a_nb   (nb = E, W, N, S)
%       M(n,n)  =  a_P^0 = rho*m.vol(n)/p.dt
%
%   t is only needed for the velocity when p.f > 0.
%   Rows that apply_bc_matrix replaces by a boundary equation (e.g. the inlet)
%   can stay empty here.
%
%   Useful pieces:
%       ux = velocity_profile(m.y, t, p);      % Ny x 1, one value per row i
%       n  = m.inner;  i = m.I(n);  j = m.J(n); % DOFs and their (i,j)
%       neighbours: E = n+1, W = n-1, N = n+m.Nx, S = n-m.Nx
%       build with triplets:  K = sparse(rows, cols, vals, m.N, m.N);

% TODO (Q7):  inner CVs, upwind convection + central-difference diffusion
% TODO (Q10): wall CVs, if you treat them as half CVs with a zero-flux face
% TODO (Q21): p.scheme == 'QUICK'
% TODO (Q22): Robin walls (p.C > 0)

error('venom_solver:todo', 'TODO: assemble K and M (Q7).');
end


function A = apply_bc_matrix(A, m, p) %#ok<INUSD>
% Put the boundary equations into the rows of A.
%
%   A comes in as M + theta*K. For every boundary DOF whose equation is NOT a
%   finite-volume balance (e.g. Dirichlet at the inlet, a CD Neumann condition),
%   put that equation in its row. The matching right-hand side values are set
%   in apply_bc_rhs.
%
%   Boundary DOF lists: m.inlet, m.outlet, m.bottom, m.top
%   (the corners are in two lists: decide which equation they get, report Q10).
%
%   Tip: overwriting rows of a sparse matrix (A(rows,:) = 0) is slow for large N.
%   If those rows are left empty in assemble_operator, simply add them:
%       A = A + sparse(rows, cols, vals, m.N, m.N);

% TODO (Q9):  inlet, Dirichlet  phi = phi_inlet(t)
% TODO (Q8):  outlet, d(phi)/dx = 0 discretized with CD
% TODO (Q10): walls,  d(phi)/dn = 0
% TODO (Q22): walls,  Robin  Gam*d(phi)/dn = C*phi

error('venom_solver:todo', 'TODO: boundary rows of A (Q8-Q10).');
end


function b = apply_bc_rhs(b, m, p, t) %#ok<INUSD>
% Right-hand side entries of the boundary equations at time t = t^{k+1}.
%
%   b comes in as (M - (1-theta)*K)*phi^k. For each row set in apply_bc_matrix,
%   overwrite b with the right-hand side of that equation.
%       phi_in = inlet_concentration(t, p);

% TODO (Q9):      inlet rows   -> phi_in
% TODO (Q8, Q10): Neumann rows -> 0

error('venom_solver:todo', 'TODO: boundary entries of b (Q8-Q10).');
end


function [phiOut, phiTot] = compute_diagnostics(phi, m, p) %#ok<INUSD>
% Scalar outputs recorded at every time step (Q14).
%   phiOut  y-averaged concentration at the outlet  (1/H) * int_0^H phi(L,y) dy  [g/L]
%   phiTot  total venom in the domain  int phi dV, per unit depth  [g/L * m^2]
%   Useful: m.outlet (DOFs of the outlet column, bottom to top), m.y, m.vol

% TODO (Q14)
phiOut = NaN;
phiTot = NaN;
end


%% =========================================================================
%  Storage
%  =========================================================================

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
