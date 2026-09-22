# Prototypes (archived, not maintained)

Coursework and exploratory scripts written before the production solver, kept
for provenance. They are **not** part of the `nucth` package, are not tested,
and use `scipy.integrate.simps`, which no longer exists in scipy >= 1.14.

The progression they document:

| File | Adds |
|---|---|
| `Schrodinger.py` | Non-relativistic radial Schrodinger equation, isotropic harmonic oscillator |
| `iho1.py` | Same, refined |
| `WS_no_spin.py` | Woods-Saxon potential, no spin-orbit |
| `WS_spin.py` | Woods-Saxon with a spin-orbit term |
| `scalar_dirac.py` | Dirac equation, scalar potential only |
| `vector_dirac.py` | Dirac equation, vector potential only |
| `dirac.py` | Both potentials, single nucleon species |

The production version of the last step is `src/nucth/dirac.py`.

`../Finite_Sphere_Quantum_Well.nb` is the Mathematica notebook that solved the
finite spherical well analytically — the closed-form cross-check the numerical
solver was first validated against.
