# Generated results

Output of the pipeline, kept so every figure in `notebooks/` and `thesis/`
reproduces without re-running the Linux-only RMF stage.

| Directory | Produced by | Columns |
|---|---|---|
| `potentials/raw_*.dat` | `external/rmf/main` (stage 1) | `r [fm], Phi, W, B, A [MeV]` |
| `potentials/*_potential.dat` | `nucth convert` | `r [fm], S_n, V_n, S_p, V_p [MeV]` |
| `energies/*_energy.txt` | `nucth structure` | `kappa, E-M [MeV]` per species |
| `densities/*_nucleon_density.dat` | `nucth structure` | `r [fm], rho_s(n), rho_v(n), rho_s(p), rho_v(p) [fm^-3]` |
| `densities/*_charge_densities.dat` | `nucth formfactors` | `r [fm], charge density, weak density [fm^-3]` |
| `densities/*_WeakFF.dat` | `nucth formfactors` | `q [fm^-1], F_weak, F_charge` |
| `cross_sections/*_cross.dat` | `nucth cevns` | `E_nu [MeV], sigma [cm^2]` |
| `skins/*_neutron_skins.txt` | `nucth formfactors` | `rp, rn, rn-rp, rch, rw, rw-rch [fm]` per model |

## Scope

Seven nuclei are committed, ~25 MB total: **Ar40** (the CEvNS target),
**Ca40**, **Ca48** and **Pb208** (the primary neutron-skin correlation), plus
**Cs133**, **I127** and **Xe132**, which `notebooks/thesis_plots.ipynb` needs
for the headline skin-vs-cross-section fit across multiple nuclei — a
three-point fit (Ar40/Ca40/Ca48/Pb208 alone) isn't the thesis result. A small
slice of **Ti50** (four RMF models' nucleon densities only) supports the
notebook's separate mirror-nuclei section.

Ni50 and Sn132 from the original working tree are not committed; they are
regenerable with

    nucth all --nucleus Sn132 --model fsugold --neutrons 82 --protons 50

given a raw potential from the RMF stage.
