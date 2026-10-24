# ME-474 Assessment 1: venom transport in a blood vessel

Group 12. Unsteady 2D convection-diffusion of a passive scalar in a plane channel,
finite volume method, theta scheme in time.

## Layout

```
default_params.m     all physical + numerical parameters (baseline = Q12)
venom_solver.m       solver + all numerics as local functions:
                       time loop  (M + theta K) phi^{k+1} = (M - (1-theta) K) phi^k
                       mesh, velocity, inlet, UD/QUICK assembly, boundary rows,
                       outlet average / total venom / leave time
plot_results.m       'timeseries' | 'snapshots' | 'field' | 'video' | 'save'
q12_q14_baseline.m   Q11-Q14
q15_q20_explicit.m   Q15, Q20
q16_q18_time_schemes.m  Q16-Q18
q19_grid.m           Q19
q21_quick.m          Q21
output/              generated figures and videos (git-ignored)
```

## Conventions (match the report)

- Nodes sit on all four boundaries: `dx = L/(Nx-1)`, `dy = H/(Ny-1)`; boundary CVs are half CVs.
- `i` = row (y), `i = 1` bottom wall, `i = Ny` top wall.
  `j` = column (x), `j = 1` inlet, `j = Nx` outlet.
- DOF numbering row by row: `n = j + (i-1)*Nx`.
  Neighbours: E = n+1, W = n-1, N = n+Nx, S = n-Nx.
- `out.mesh.toField(phi)` turns a DOF vector into an `Ny x Nx` array `Phi(i,j)`.

## Workflow

Run everything from the repo root.

Each script starts from `default_params()` and changes only what the question asks.

## Progress

| Questions | Topic | Status |
|---|---|---|
| Q7 | inner CV assembly | done |
| Q8-Q10 | outlet, inlet, wall BCs | done |
| Q12-Q14 | baseline run, BC checks, outlet/total venom | done |
| Q15-Q18 | time step and theta studies | done |
| Q19-Q20 | grid refinement | done |
| Q21 | QUICK | done |
| Q22 | Robin (absorbing) walls | todo |
| Q23-Q24 | pulsating flow, Pe = 100, video | todo |
