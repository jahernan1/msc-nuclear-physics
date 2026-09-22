# msc-nuclear-physics

**Coherent elastic neutrino-nucleus scattering (CEvNS) cross sections from
relativistic mean-field nuclear structure.**

MSc thesis code, Florida State University, 2019. The physics question: how
sensitive is the CEvNS cross section to the **neutron skin** — the difference
between a nucleus's neutron and proton radii — across relativistic mean-field
models that are otherwise equally good at reproducing known nuclear data?

**Thesis:** [Jesse A. Hernandez, MSc thesis, Florida State University, 2019](https://purl.lib.fsu.edu/diginole/2019_Fall_Hernandez_fsu_0071N_15473)
**Published paper:** Junjie Yang, Jesse A. Hernandez, and J. Piekarewicz, ["Electroweak probes of ground state densities"](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.100.054301), Phys. Rev. C **100**, 054301 (2019)

## The pipeline

```

   RMF parameter set  (fsugold, fsugarnet, rmf022, tamufsua, ...)
            |
   [1] RMF mean-field solve                    external/rmf/main  (Linux binary)
            |   sigma, omega, rho, photon fields
            v
       raw_*.dat                   r, Phi, W, B, A  [MeV]
            |
   [2] meson fields -> nucleon potentials      nucth convert
            |   S(p)=S(n)=-Phi   V(n)=W-B/2   V(p)=W+B/2+A
            v*
       _potential.dat
            |
   [3] coupled radial Dirac equations          nucth structure
            |   shooting + matching, RK4, Brent root-finding
            v
       single-particle spectrum + scalar/vector densities
            |
   [4] Fourier transform, fold nucleon FF      nucth formfactors
            |   F_ch = Gep F_p + Gen F_n
            |   F_w  = Gwp F_p + Gwn F_n
            v
       charge & weak densities, r_p, r_n, r_ch, r_w, neutron skin
            |
   [5] CEvNS cross section                     nucth cevns
            |   dsigma/dT = (Gf^2/4pi) M Qw^2 [1 - T/E - MT/2E^2] |F_w(Q)|^2
            v
       sigma(E_nu)  [cm^2]

```

**Why the neutron skin matters here:** the weak charge is
`Qw = N - (1 - 4 sin^2 theta_W) Z`, and `1 - 4 sin^2 theta_W` is approximately
0.076. The proton contribution is suppressed by more than an order of
magnitude, so CEvNS is a *neutron*-distribution probe. The weak form factor
`F_w(Q)` carries that dependence, and the cross section inherits it.

## The physics core

The coupled radial Dirac equations for the upper and lower components `g, f` of
a nucleon in scalar potential `S` and vector potential `V`:

```

dg/dx = (E - V + S + 1) f - (kappa/x) g
df/dx = (kappa/x) f - (E - V - S - 1) g

```

in dimensionless form (`x = r M / hbarc`, energies in units of `M`). `kappa` is
the Dirac quantum number: `j = |kappa| - 1/2`, and `l = kappa` for `kappa > 0`,
`l = -1 - kappa` for `kappa < 0`.

Bound states come from shooting and matching: RK4 outward from the origin,
inward from the box boundary, then Brent root-finding on the determinant
`g_out f_in - g_in f_out`. Levels fill to `N` and `Z` by ascending energy, giving

```

rho_v(r) = sum_occ (2j+1) / (4 pi r^2) * (g^2 + f^2)
rho_s(r) = sum_occ (2j+1) / (4 pi r^2) * (g^2 - f^2)

```

## Install

With [uv](https://docs.astral.sh/uv/) (recommended):

```bash
git clone https://github.com/<you>/msc-nuclear-physics
cd msc-nuclear-physics
uv venv --python 3.11
source .venv/bin/activate
uv pip install -e ".[dev,plots]"
```

With conda:

```bash
conda env create -f environment.yml
conda activate nucth
```

## Run

Stages 2 through 5, for Ca40 in the FSUGold parameterization:

```bash
nucth all --nucleus Ca40 --model fsugold --neutrons 20 --protons 20
```

Or one stage at a time:

```bash
nucth convert     --nucleus Ca40 --model fsugold
nucth structure   --nucleus Ca40 --model fsugold --neutrons 20 --protons 20
nucth formfactors --nucleus Ca40 --model fsugold
nucth cevns       --nucleus Ca40 --model fsugold --neutrons 20 --protons 20
```

Heavy nuclei need a larger matching radius and box:

```bash
nucth structure --nucleus Pb208 --model fsugold --neutrons 126 --protons 82 \
                --r-match 10 --r-box 20 --kappa-max 8
```

## Reproducibility boundary

**Stage 1 is not reproducible from this repository.** The RMF solver survives
only as a compiled Linux x86-64 binary — no source was preserved. See
`external/rmf/README.md`.

**Stages 2 through 5 are fully reproducible** from the committed
`results/potentials/*_potential.dat` files on any platform, and the test suite
pins them to the numbers that went into the thesis.

## Tests

```bash
pytest -m "not slow"   # unit tests, seconds
pytest -m slow         # characterization tests, minutes
pytest                 # everything
```

The slow tests are the interesting ones: they run the full Dirac solver on
`results/potentials/Ca40_fsugold_potential.dat` and assert that the
single-particle spectrum matches `results/energies/Ca40_fsugold_energy.txt` and
that the radii and skins match `results/skins/Ca40_neutron_skins.txt` — both
files produced by the original 2019 code. They are what makes the refactor
verifiable rather than merely plausible.

## Layout

| Path                | Contents                                                                        |
| ------------------- | ------------------------------------------------------------------------------- |
| `src/nucth/`        | The package: `transforms`, `potentials`, `dirac`, `formfactors`, `cevns`, `cli` |
| `tests/`            | Unit and characterization tests                                                 |
| `data/`             | Nucleon form factors (`GPN.txt`), nuclear masses, 13 RMF parameter sets         |
| `external/rmf/`     | The stage-1 binary and its documentation                                        |
| `results/`          | Generated output for Ar40, Ca40, Ca48, Cs133, I127, Pb208, Xe132 (+ a Ti50 slice) |
| `notebooks/`        | Thesis and form-factor figure notebooks                                         |
| `figures/`          | Final thesis figures                                                            |
| `thesis/`           | LaTeX source and the compiled thesis                                            |
| `archive/`          | Pre-production prototypes, kept for provenance                                  |

## References

The physics background is in `thesis/myrefs.bib`. Key external results this
work builds on: the COHERENT collaboration's first CEvNS observation (2017),
PREX parity-violating electron scattering on Pb208, and Horowitz & Piekarewicz
on the neutron-skin/symmetry-energy correlation.

## License

MIT for the Python code — see `LICENSE`. `external/rmf/main` is a third-party
compiled binary and is not covered by it.
