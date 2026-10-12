# ME-474 Assessment 1: venom transport in a blood vessel

Group 12. Unsteady 2D convection-diffusion of a passive scalar in a plane channel,
finite volume method, theta scheme in time.

## Layout

```
default_params.m     all physical + numerical parameters (baseline = Q12)
venom_solver.m       solver + all numerics as local functions:
                       time loop  (M + theta K) phi^{k+1} = (M - (1-theta) K) phi^k
                       build_mesh, velocity_profile, inlet_concentration
                       TODO: assemble_operator (Q7, Q10, Q21, Q22),
                             apply_bc_matrix / apply_bc_rhs (Q8-Q10, Q22),
                             compute_diagnostics (Q14)
plot_results.m       'timeseries' | 'snapshots' | 'field' | 'video' | 'save'
qXX_*.m              one run script per question group, e.g. q12_q14_baseline.m
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

1. Fill in the TODO local functions in `venom_solver.m`.
2. Run `q12_q14_baseline.m`, then add one script per question
   (copy `default_params`, change only what the question asks).

## Progress

| Questions | Topic | Status |
|---|---|---|
| Q7 | inner CV assembly | todo |
| Q8-Q10 | outlet, inlet, wall BCs | todo |
| Q12-Q14 | baseline run, BC checks, outlet/total venom | todo |
| Q15-Q18 | time step and theta studies | todo |
| Q19-Q20 | grid refinement | todo |
| Q21 | QUICK | todo |
| Q22 | Robin (absorbing) walls | todo |
| Q23-Q24 | pulsating flow, Pe = 100, video | todo |
