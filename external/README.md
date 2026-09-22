# RMF mean-field stage (external, binary only)

Stage 1 of the pipeline — the self-consistent relativistic mean-field solve that
produces the meson fields (sigma, omega, rho, photon) for a given nucleus and
parameter set — was run with the `main` binary in this directory.

**No source code for** `main` **survives in this repository.** It was a C/C++ program
built against GSL 2.3 by the FSU nuclear theory group; only the compiled artifact
and its build notes (`insturction.txt`) were preserved for this project.

`main` is a **Linux x86-64 ELF executable**. It does not run on macOS or Windows
without a container or VM.

## Invocation

```
./main <parameter_set> <Z> <N> <label> <kappa_max> [<n_extra>]
```

Example, from `script.sh`:

```
./main tamufsua 18 22 othermodels 5
```

Outputs, written next to the binary:

- `potential0.txt` — columns: `r [fm], Phi [MeV], W [MeV], B [MeV], A [MeV]`
- `density0.txt`   — columns: `r [fm], rho_s(n), rho_v(n), rho_s(p), rho_v(p)`
- `output/` — total energies, single-state energies, radii

`potential0.txt` is what enters this repository as
`results/potentials/raw_<nucleus>_<model>.dat`, and `nucth.potentials.convert`
turns it into the nucleon scalar/vector potentials the Dirac solver consumes.

## Reproducibility boundary

**Everything downstream of stage 1 is fully reproducible from the committed**
`results/potentials/*_potential.dat` **files using this Python package.** Stage 1
itself is reproducible only on Linux x86-64 with this binary.

## Parameter sets

13 sets are committed under `data/parameters/`: fsu020, fsugold2, iufsu,
parameter_garnet, parameter_gold, parameter_linear, parameter_nl3, rmf022,
rmf028, rmf032, tamufsua, tamufsub, tamufsuc.

**Known gap:** potentials for `rmf012` and `rmf016` are committed under
`results/potentials/`, but their parameter files were not preserved. Those two
models cannot be regenerated from scratch.

## Provenance of `GPN.txt`

`NuclearFF.f` is the Fortran program that generated `data/GPN.txt`
(columns: `q [fm^-1], Gep, Gen` — 8001 rows), the single-nucleon electric form
factors folded into the nuclear form factor in `nucth.formfactors`.