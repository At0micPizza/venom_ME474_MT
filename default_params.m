function p = default_params()
% DEFAULT_PARAMS  Baseline setup of the venom problem (Q11-Q12).
%   Copy the struct in a run script and change only what the question asks:
%       p = default_params();  p.theta = 1;  out = venom_solver(p);

%% physics (Q11)
p.L     = 1;          % channel length [m]
p.H     = 0.01;       % channel height [m]
p.rho   = 1050;       % blood density [kg/m^3]
p.nu    = 3.5e-6;     % kinematic viscosity [m^2/s]
p.umean = 0.1;        % mean streamwise velocity [m/s]
p.Pe    = 10;         % Peclet number  Pe = rho*umean*H/Gam  -> sets Gam
p.Gam   = [];         % diffusion coeff [kg/(m s)]; empty = computed from Pe

%% venom injection (Q12)
p.phi0    = 0;        % initial concentration in the whole domain [g/L]
p.tBite   = 0.1;      % bite duration [s]
p.phiBite = 1;        % inlet concentration during the bite [g/L]

%% walls (Q22)
p.C = 0;              % Robin absorption coefficient [kg/(m^2 s)]; 0 = non-absorbing

%% pulsating flow (Q23-Q24)
p.f     = 0;          % heart frequency [Hz]; 0 = steady Poiseuille flow
p.alpha = [];         % amplitude factor of the pulsating flow (Q23); set it in the run script

%% numerics
p.Nx     = 100;       % CVs in x (nodes on inlet and outlet)
p.Ny     = 10;        % CVs in y (nodes on both walls)
p.dt     = 0.005;     % time step [s]
p.tEnd   = 18;        % end time [s]
p.theta  = 0;         % 0 = explicit Euler, 1 = implicit Euler, 0.5 = Crank-Nicolson
p.scheme = 'UD';      % convection scheme: 'UD' or 'QUICK' (Q21)

%% output
p.name          = 'baseline';  % label used in legends and file names
p.snapshotTimes = [];          % store the full 2D field at these times [s] (Q13)
p.frameEvery    = 0;           % store a video frame every n steps; 0 = no frames
p.divergeLimit  = 1e3;         % stop the run when max|phi| exceeds this
p.verbose       = true;
end
