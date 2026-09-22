"""Command-line pipeline.

    nucth convert     --nucleus Ca40 --model fsugold
    nucth structure   --nucleus Ca40 --model fsugold --neutrons 20 --protons 20
    nucth formfactors --nucleus Ca40 --model fsugold
    nucth cevns       --nucleus Ca40 --model fsugold --neutrons 20 --protons 20
    nucth all         --nucleus Ca40 --model fsugold --neutrons 20 --protons 20

This is the only module that prints, plots, or writes files as a side effect.
Replaces `nuc-th/allNuclear.py`, which hardcoded Dropbox paths and re-invoked
each stage through `ipython` as a subprocess.

Stage 1 (the RMF solve itself) is not driven from here: it is a Linux-only
binary, documented in `external/rmf/README.md`. Start from `convert` with a
raw potential already in `results/potentials/`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from nucth import cevns as cevns_mod
from nucth import formfactors as ff_mod
from nucth import potentials as pot_mod
from nucth.dirac import Nuclear, binding_energies

DEFAULT_RESULTS = Path("results")

NUCLEON_DENSITY_HEADER = "raxis[fm], rho_s(n), rho_v(n), rho_s(p), rho_v(p) [fm^-3]"
WEAK_FF_HEADER = "q[fm^-1], F_weak, F_charge [normalized to 1 at q=0]"
CHARGE_DENSITY_HEADER = "r[fm], charge density, weak density [fm^-3]"
CROSS_HEADER = "E_nu [MeV], sigma [cm^2]"


def run_convert(raw_path: Path, out_path: Path) -> None:
    """Stage 2: raw meson fields -> nucleon scalar/vector potentials."""
    potential = pot_mod.convert(pot_mod.load_raw(raw_path))
    pot_mod.save_potential(out_path, potential)
    print(f"wrote {out_path}")


def run_structure(
    potential_path: Path,
    N: int,
    Z: int,
    kappa_max: int,
    energy_out: Path,
    density_out: Path,
    *,
    r_match: float = 7.0,
    r_box: float = 15.0,
) -> None:
    """Stage 3: solve the Dirac equations, write the spectrum and densities."""
    nucleus = Nuclear(potential_path, N=N, Z=Z, r_match=r_match, r_box=r_box)
    kappa_values = list(range(-kappa_max, 0)) + list(range(1, kappa_max))

    neutron_levels, proton_levels = nucleus.spectrum(kappa_values)
    n_occ = nucleus.occupied("neutron", neutron_levels)
    p_occ = nucleus.occupied("proton", proton_levels)

    energy_out.parent.mkdir(parents=True, exist_ok=True)
    with energy_out.open("w") as handle:
        handle.write(f"#{'kappa':_<5}{'Neutron':_^25}{'kappa':_<5}{'Proton':_^25}\n")
        n_mev = binding_energies(n_occ, nucleus.MASS)
        p_mev = binding_energies(p_occ, nucleus.MASS)
        for i in range(max(len(n_occ), len(p_occ))):
            nk, ne = (n_occ[i][0], n_mev[i]) if i < len(n_occ) else (0, 0.0)
            pk, pe = (p_occ[i][0], p_mev[i]) if i < len(p_occ) else (0, 0.0)
            handle.write(f"{nk:< 5}{ne:^25}{pk:< 5}{pe:^25}\n")
    print(f"wrote {energy_out}")

    n_s, n_v = nucleus.density("neutron", nucleus.wavefunctions("neutron", n_occ))
    p_s, p_v = nucleus.density("proton", nucleus.wavefunctions("proton", p_occ))

    density_out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        density_out,
        np.stack((nucleus.radial_axis, n_s, n_v, p_s, p_v), axis=-1),
        header=NUCLEON_DENSITY_HEADER,
    )
    print(f"wrote {density_out}")


def run_formfactors(density_path: Path, charge_out: Path, weakff_out: Path) -> None:
    """Stage 4: charge and weak form factors, densities, radii, skins."""
    result = ff_mod.compute(density_path)

    print(f"electric charge  Q_ch = {result.Q_charge:.5f}")
    print(f"weak charge      Q_w  = {result.Q_weak:.5f}")
    print(f"proton radius    r_p  = {result.r_p:.5f} fm")
    print(f"neutron radius   r_n  = {result.r_n:.5f} fm")
    print(f"neutron skin           = {result.neutron_skin:+.5f} fm")
    print(f"charge radius    r_ch = {result.r_ch:.5f} fm")
    print(f"weak radius      r_w  = {result.r_w:.5f} fm")
    print(f"weak skin              = {result.weak_skin:+.5f} fm")

    charge_out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        charge_out,
        np.stack((result.r, result.charge_density, result.weak_density), axis=-1),
        header=CHARGE_DENSITY_HEADER,
    )
    np.savetxt(
        weakff_out,
        np.stack((result.q, result.weak_ff, result.charge_ff), axis=-1),
        header=WEAK_FF_HEADER,
    )
    print(f"wrote {charge_out}\nwrote {weakff_out}")


def run_cevns(
    weakff_path: Path,
    N: int,
    Z: int,
    cross_out: Path,
    *,
    e_max: float = 52.8,
    n_points: int = 500,
) -> None:
    """Stage 5: CEvNS cross section as a function of neutrino energy."""
    form_factor = cevns_mod.WeakFormFactor.from_file(weakff_path)
    mass = cevns_mod.nuclear_mass(N, Z)
    qw = cevns_mod.weak_charge(N, Z)
    print(f"nuclear mass M = {mass:.3f} MeV, weak charge Qw = {qw:.4f}")

    t_min = 1.0e-6
    e_min = 0.5 * (t_min + np.sqrt(t_min**2 + 2.0 * mass * t_min)) + 0.01
    energies = np.linspace(e_min, e_max, n_points)
    sigma = np.array(
        [cevns_mod.total_cross_section(mass, e, qw, form_factor) for e in energies]
    )

    cross_out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(cross_out, np.stack((energies, sigma), axis=-1), header=CROSS_HEADER)
    print(f"sigma({energies[-1]:.1f} MeV) = {sigma[-1]:.4e} cm^2")
    print(f"wrote {cross_out}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="nucth", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    for name in ("convert", "structure", "formfactors", "cevns", "all"):
        p = sub.add_parser(name)
        p.add_argument("--nucleus", required=True, help="e.g. Ca40")
        p.add_argument("--model", required=True, help="e.g. fsugold")
        p.add_argument("--neutrons", type=int, default=0)
        p.add_argument("--protons", type=int, default=0)
        p.add_argument("--kappa-max", type=int, default=5)
        p.add_argument("--r-match", type=float, default=7.0, help="matching radius [fm]")
        p.add_argument("--r-box", type=float, default=15.0, help="box boundary [fm]")
        p.add_argument("--results", type=Path, default=DEFAULT_RESULTS)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    tag = f"{args.nucleus}_{args.model}"
    root = args.results

    raw = root / "potentials" / f"raw_{tag}.dat"
    potential = root / "potentials" / f"{tag}_potential.dat"
    energy = root / "energies" / f"{tag}_energy.txt"
    density = root / "densities" / f"{tag}_nucleon_density.dat"
    charge = root / "densities" / f"{tag}_charge_densities.dat"
    weakff = root / "densities" / f"{tag}_WeakFF.dat"
    cross = root / "cross_sections" / f"{tag}_cross.dat"

    stages = [args.command] if args.command != "all" else [
        "convert", "structure", "formfactors", "cevns"
    ]

    for stage in stages:
        if stage == "convert":
            run_convert(raw, potential)
        elif stage == "structure":
            run_structure(
                potential, args.neutrons, args.protons, args.kappa_max,
                energy, density, r_match=args.r_match, r_box=args.r_box,
            )
        elif stage == "formfactors":
            run_formfactors(density, charge, weakff)
        elif stage == "cevns":
            run_cevns(weakff, args.neutrons, args.protons, cross)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
