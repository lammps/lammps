#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the *s-band* overlay file for a two-band ("double-band") EAM
potential, for use with::

    pair_style hybrid/overlay eam/fs eam/fs
    pair_coeff * * eam/fs 1  <dband>.eam.fs  Fe Cr   # pair V + d-band  F_d(rho_d)
    pair_coeff * * eam/fs 2  <sband>.eam.fs  Fe Cr   # s-band  F_s(rho_s), no pair

Background
----------
The two-band model (Olsson et al., Phys. Rev. B 72, 214119 (2005); Bonny et
al., Phil. Mag. 91, 1724 (2011)) writes the energy of atom i as

    E_i = 1/2 sum_j V(r_ij) + F_d(rho_d,i) + F_s(rho_s,i)

i.e. an ordinary EAM with an *extra* embedding term F_s acting on a second
("s-band") density.  Grouping the terms,

    E_i = [ 1/2 sum_j V + F_d(rho_d) ]  +  [ F_s(rho_s) ]
            \\___ a standard eam/fs ___/      \\__ this file: eam/fs, V == 0 __/

makes the whole model an exact `hybrid/overlay` of two `eam/fs` potentials,
needing no change to LAMMPS.  See doc/src/Howto_eam_overlay.rst.

The s-band density vanishes for *like* pairs and is non-zero only for
*unlike* pairs (Bonny Eq. 3: phi_s_AA = phi_s_BB = 0, phi_s_AB = phi_s_BA).
Such a cross-only density is expressible only in the Finnis-Sinclair
(`eam/fs`) form, where the density is indexed per ordered element pair --
it cannot be written as `eam/alloy`.  That is why both overlay files use
`eam/fs` and why this script emits an `eam/fs` file.

What this script does (and does not) guarantee
----------------------------------------------
The file *format* and the cross-only/zero-pair *structure* are exact and are
the genuinely reusable part of this tool.  The numeric s-band parameters that
reproduce a *published* Fe-Cr potential (the absolute normalisation of rho_s
and the matching F_s coefficients) are subtle and convention-dependent; the
PRESETS below collect the literature values as documented starting points but
they should be validated against the source (e.g. the mixing-enthalpy curve)
before production use.  The default `--demo` parameters are self-consistent
and *illustrative* only: they correctly exercise the overlay mechanics and
isolate the concentration-dependent s-band term, but are not a published fit.

Usage examples
--------------
    # Stand-alone, runnable two-band demo (writes BOTH files + a sample input):
    ./two_band_fecr.py --demo --prefix FeCr_demo

    # Real use: add an s-band to your own trusted Fe-Cr d-band eam/fs file:
    ./two_band_fecr.py --dband MyFeCr.eam.fs --preset olsson-vasp --prefix FeCr_2bm
"""

import argparse
import math
import os
import sys
from datetime import date

# --------------------------------------------------------------------------
# s-band parameter presets.
#
# slater: phi_s(r) = K * r**6 * exp(-2*zeta*r) * switch(r; rc_in, rc_out)
#         (the square of an r^3 4s-Slater orbital; Olsson Eq. 7 / Bonny Eq. 8).
#         phi_s is applied to UNLIKE pairs only; like pairs are identically 0.
# fs:     F_s(rho) = c1*sqrt(rho) + c2*rho**2 + c3*rho**4   (Olsson Eq. 2 form)
#         given per element (here equal for both, as in Olsson).
#
# NOTE: see the docstring -- the absolute rho_s normalisation (K) and the
# matching F_s coefficients are convention-sensitive.  Treat as starting points.
# --------------------------------------------------------------------------
PRESETS = {
    # Olsson 2005, fit to VASP mixing enthalpy (their Table III, "VASP").
    "olsson-vasp": {
        "slater": dict(K=25.0, zeta=1.323, rc_in=5.0, rc_out=5.3),
        "fs": dict(c1=-0.503, c2=-0.60, c3=0.50),
        "note": "Olsson PRB 72 214119 (2005), Table III (VASP fit); Eq. 7 Slater.",
    },
    # Olsson 2005, fit to EMTO mixing enthalpy (their Table III, "EMTO").
    "olsson-emto": {
        "slater": dict(K=25.0, zeta=1.323, rc_in=5.0, rc_out=5.3),
        "fs": dict(c1=-0.800, c2=-1.00, c3=0.80),
        "note": "Olsson PRB 72 214119 (2005), Table III (EMTO fit); Eq. 7 Slater.",
    },
    # Bonny 2011 Slater (Eq. 8); F_s = c1*sqrt(rho)+c2*rho**2 (Eq. 10).
    # The c1,c2 below are PLACEHOLDERS -- supply from Bonny et al. Table.
    "bonny": {
        "slater": dict(K=20.34075425, zeta=2.5001, rc_in=5.1, rc_out=5.3),
        "fs": dict(c1=-1.0, c2=0.0, c3=0.0),
        "note": "Bonny PhilMag 91 1724 (2011), Eq. 8 Slater; F_s coeffs are "
                "PLACEHOLDERS -- replace c1,c2 with the published values.",
    },
    # Self-consistent illustrative parameters (NOT a published fit).  Chosen so
    # rho_s is O(1) and F_s contributes O(0.05 eV) at alloyed sites.
    "demo": {
        "slater": dict(K=40.0, zeta=2.5, rc_in=5.0, rc_out=5.3),
        "fs": dict(c1=-0.30, c2=0.40, c3=0.0),
        "note": "Illustrative only -- exercises the overlay, not a real fit.",
    },
}

FLOATS_PER_LINE = 5


# --------------------------------------------------------------------------
# analytic s-band functions
# --------------------------------------------------------------------------
def switch(r, rc_in, rc_out):
    """Smooth sine switching: 1 for r<=rc_in, 0 for r>=rc_out (Bonny Eq. 9)."""
    if r <= rc_in:
        return 1.0
    if r >= rc_out:
        return 0.0
    rm = 0.5 * (rc_in + rc_out)
    d = 0.5 * (rc_out - rc_in)
    return 0.5 * (1.0 - math.sin(0.5 * math.pi * (r - rm) / d))


def slater_density(r, K, zeta, rc_in, rc_out):
    """Cross s-band density phi_s(r) = K r^6 exp(-2 zeta r) * switch(r)."""
    if r <= 0.0:
        return 0.0
    return K * r**6 * math.exp(-2.0 * zeta * r) * switch(r, rc_in, rc_out)


def fs_embedding(rho, c1, c2, c3):
    """F_s(rho) = c1*sqrt(rho) + c2*rho^2 + c3*rho^4."""
    if rho < 0.0:
        rho = 0.0
    return c1 * math.sqrt(rho) + c2 * rho * rho + c3 * rho**4


# --------------------------------------------------------------------------
# eam/fs writer
# --------------------------------------------------------------------------
def fmt_block(values):
    """Format a list of floats as text lines, FLOATS_PER_LINE per line."""
    out = []
    for i in range(0, len(values), FLOATS_PER_LINE):
        chunk = values[i:i + FLOATS_PER_LINE]
        out.append(" " + " ".join("% .16E" % v for v in chunk))
    return "\n".join(out)


def write_eam_fs(path, comments, elements, masses, znums, a0s, lats,
                 nrho, drho, nr, dr, cut, frho, rhor, z2r):
    """Write a DYNAMO setfl Finnis-Sinclair (eam/fs) file.

    frho[i][m]    : embedding F_i(rho), m = 0..nrho-1 (rho = m*drho)
    rhor[i][j][m] : density emitted by element i, received by element j
                    (LAMMPS indexes rho at receiver j from neighbor i)
    z2r[i][j][m]  : r*phi_ij(r), lower triangle i>=j
    """
    nel = len(elements)
    with open(path, "w") as fp:
        # three comment lines (line 1 follows the LAMMPS metadata convention)
        for line in comments[:3]:
            fp.write(line.rstrip() + "\n")
        fp.write("%d %s\n" % (nel, " ".join(elements)))
        fp.write("%d % .16E %d % .16E % .16E\n" % (nrho, drho, nr, dr, cut))
        for i in range(nel):
            fp.write("%d % .16E % .16E %s\n" % (znums[i], masses[i], a0s[i], lats[i]))
            fp.write(fmt_block(frho[i]) + "\n")
            for j in range(nel):
                fp.write(fmt_block(rhor[i][j]) + "\n")
        for i in range(nel):
            for j in range(i + 1):
                fp.write(fmt_block(z2r[i][j]) + "\n")


# --------------------------------------------------------------------------
# minimal eam/fs reader (used by --demo to reuse a shipped 1-element Fe block)
# --------------------------------------------------------------------------
def _read_floats(lines, pos, count):
    """Consume whole lines from `lines` starting at `pos` until `count`
    floats are collected.  Returns (values, new_pos)."""
    vals = []
    while len(vals) < count:
        vals.extend(float(t) for t in lines[pos].split())
        pos += 1
    if len(vals) != count:
        raise ValueError("misaligned numeric block (expected %d, got %d)"
                         % (count, len(vals)))
    return vals, pos


def read_eam_fs(path):
    """Read an eam/fs file into a dict.  Handles the general Nelement format."""
    with open(path) as fp:
        lines = [ln for ln in fp.read().splitlines()]
    tok = lines[3].split()
    nel = int(tok[0])
    elements = tok[1:1 + nel]
    g = lines[4].split()
    nrho, drho, nr, dr, cut = int(g[0]), float(g[1]), int(g[2]), float(g[3]), float(g[4])
    masses, znums, a0s, lats = [], [], [], []
    frho = [None] * nel
    rhor = [[None] * nel for _ in range(nel)]
    pos = 5
    for i in range(nel):
        hdr = lines[pos].split()
        pos += 1
        znums.append(int(float(hdr[0])))
        masses.append(float(hdr[1]))
        a0s.append(float(hdr[2]) if len(hdr) > 2 else 0.0)
        lats.append(hdr[3] if len(hdr) > 3 else "unknown")
        frho[i], pos = _read_floats(lines, pos, nrho)
        for j in range(nel):
            rhor[i][j], pos = _read_floats(lines, pos, nr)
    z2r = [[None] * nel for _ in range(nel)]
    for i in range(nel):
        for j in range(i + 1):
            z2r[i][j], pos = _read_floats(lines, pos, nr)
    return dict(elements=elements, masses=masses, znums=znums, a0s=a0s, lats=lats,
                nrho=nrho, drho=drho, nr=nr, dr=dr, cut=cut,
                frho=frho, rhor=rhor, z2r=z2r)


# --------------------------------------------------------------------------
# builders
# --------------------------------------------------------------------------
def build_sband(elements, masses, znums, a0s, lats, pair, params,
                nrho, drho, nr, dr, cut, citation):
    """Build the s-band eam/fs arrays (cross-only density, F_s, zero pair)."""
    nel = len(elements)
    ia, ib = pair
    sp = params["slater"]
    fsc = params["fs"]

    frho = [[fs_embedding(m * drho, **fsc) for m in range(nrho)] for _ in range(nel)]

    # density emitted by i, received by j: non-zero ONLY for the unlike pair.
    rhor = [[[0.0] * nr for _ in range(nel)] for _ in range(nel)]
    cross = [slater_density(m * dr, **sp) for m in range(nr)]
    rhor[ia][ib] = list(cross)   # A emits to B
    rhor[ib][ia] = list(cross)   # B emits to A (symmetric)

    # no pair term in the s-band file
    z2r = [[[0.0] * nr for _ in range(nel)] for _ in range(nel)]

    comments = [
        "DATE: %s UNITS: metal CONTRIBUTOR: two_band_fecr.py "
        "CITATION: %s" % (date.today().isoformat(), citation),
        "two-band EAM s-band overlay file (use with pair_style hybrid/overlay "
        "eam/fs eam/fs)",
        "cross-only s-density between %s-%s; pair term is zero; %s"
        % (elements[ia], elements[ib], params["note"]),
    ]
    return dict(comments=comments, elements=elements, masses=masses, znums=znums,
                a0s=a0s, lats=lats, nrho=nrho, drho=drho, nr=nr, dr=dr, cut=cut,
                frho=frho, rhor=rhor, z2r=z2r)


def build_demo_dband(fe_file, elements):
    """Build an *illustrative* 2-element d-band eam/fs by duplicating the
    Fe block from a shipped 1-element Fe eam/fs file for both elements.

    The two elements are physically identical here, so the ONLY thing that
    distinguishes them is the s-band overlay -- which cleanly isolates and
    demonstrates the concentration-dependent term.  This is a teaching aid,
    NOT a real Fe-Cr potential.
    """
    fe = read_eam_fs(fe_file)
    if len(fe["elements"]) != 1:
        sys.exit("--demo expects a single-element Fe eam/fs file (e.g. "
                 "potentials/Fe_mm.eam.fs); got %d elements" % len(fe["elements"]))
    nr, nrho = fe["nr"], fe["nrho"]
    nel = len(elements)
    masses = [fe["masses"][0]] * nel
    znums = [fe["znums"][0]] * nel
    a0s = [fe["a0s"][0]] * nel
    lats = [fe["lats"][0]] * nel
    frho = [list(fe["frho"][0]) for _ in range(nel)]
    # d-density depends only on the emitter: every element emits the Fe density.
    rhor = [[list(fe["rhor"][0][0]) for _ in range(nel)] for _ in range(nel)]
    # pair term: same V for all pairs (lower triangle filled)
    z2r = [[list(fe["z2r"][0][0]) for _ in range(nel)] for _ in range(nel)]
    comments = [
        "DATE: %s UNITS: metal CONTRIBUTOR: two_band_fecr.py "
        "CITATION: illustrative demo derived from %s"
        % (date.today().isoformat(), os.path.basename(fe_file)),
        "ILLUSTRATIVE two-band demo d-band: both elements are physically Fe",
        "only the s-band overlay distinguishes the elements; NOT a real Fe-Cr fit",
    ]
    return dict(comments=comments, elements=elements, masses=masses, znums=znums,
                a0s=a0s, lats=lats, nrho=nrho, drho=fe["drho"], nr=nr, dr=fe["dr"],
                cut=fe["cut"], frho=frho, rhor=rhor, z2r=z2r)


def write_sample_input(path, dband, sband, elements):
    el = " ".join(elements)
    text = """\
# Two-band (s+d) EAM via hybrid/overlay -- sample input (illustrative demo).
# Pure-element regions feel only the d-band; mixed regions also feel F_s(rho_s).

units           metal
atom_style      atomic
boundary        p p p

lattice         bcc 2.8553
region          box block 0 6 0 6 0 6
create_box      2 box
create_atoms    1 box
# make ~half the atoms type 2 to create an alloy and a non-zero s-band energy
set             group all type/fraction 2 0.5 12345

mass            1 55.845
mass            2 55.845

pair_style      hybrid/overlay eam/fs eam/fs
pair_coeff      * * eam/fs 1 {dband} {el}
pair_coeff      * * eam/fs 2 {sband} {el}

# decompose the energy: total, d-band sub-style, s-band sub-style
compute         ed all pair eam/fs 1 epair
compute         es all pair eam/fs 2 epair
thermo_style    custom step temp pe c_ed c_es press
thermo          10

velocity        all create 300.0 87287 loop geom
fix             1 all nve
run             100
""".format(dband=dband, sband=sband, el=el)
    with open(path, "w") as fp:
        fp.write(text)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate the s-band overlay eam/fs file for a two-band "
                    "(double-band) EAM potential.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("--prefix", default="FeCr_2bm",
                    help="output filename prefix (default: FeCr_2bm)")
    ap.add_argument("--elements", nargs="+", default=["Fe", "Cr"],
                    help="element names in file order (default: Fe Cr)")
    ap.add_argument("--preset", default="demo", choices=sorted(PRESETS),
                    help="s-band parameter preset (default: demo)")
    ap.add_argument("--pair", nargs=2, metavar=("A", "B"),
                    help="the unlike element pair carrying the s-density "
                         "(default: first two elements)")
    ap.add_argument("--masses", nargs="+", type=float,
                    help="per-element masses (default: 55.845 51.996 ... )")
    ap.add_argument("--nrho", type=int, default=10000)
    ap.add_argument("--drho", type=float, default=1.0e-2)
    ap.add_argument("--nr", type=int, default=10000)
    ap.add_argument("--dr", type=float, default=5.3e-4)
    ap.add_argument("--cut", type=float, default=5.3)
    ap.add_argument("--dband", help="existing Fe-Cr d-band eam/fs file (real use)")
    ap.add_argument("--demo", action="store_true",
                    help="also emit an illustrative demo d-band + sample input "
                         "so the overlay runs end-to-end")
    ap.add_argument("--fe-file", default=None,
                    help="single-element Fe eam/fs file for --demo "
                         "(default: ../../potentials/Fe_mm.eam.fs)")
    args = ap.parse_args()

    elements = args.elements
    nel = len(elements)
    if nel < 2:
        ap.error("need at least two elements for a two-band potential")

    default_mass = {"Fe": 55.845, "Cr": 51.996, "Ni": 58.693, "Cu": 63.546}
    masses = args.masses or [default_mass.get(e, 1.0) for e in elements]
    default_z = {"Fe": 26, "Cr": 24, "Ni": 28, "Cu": 29}
    znums = [default_z.get(e, 1) for e in elements]
    a0s = [2.8553] * nel
    lats = ["bcc"] * nel

    if args.pair:
        try:
            pair = (elements.index(args.pair[0]), elements.index(args.pair[1]))
        except ValueError:
            ap.error("--pair elements must be among --elements")
    else:
        pair = (0, 1)

    params = PRESETS[args.preset]
    citation = params["note"]

    sband = build_sband(elements, masses, znums, a0s, lats, pair, params,
                        args.nrho, args.drho, args.nr, args.dr, args.cut, citation)
    sband_path = args.prefix + ".sband.eam.fs"
    write_eam_fs(sband_path, sband["comments"], sband["elements"], sband["masses"],
                 sband["znums"], sband["a0s"], sband["lats"], sband["nrho"],
                 sband["drho"], sband["nr"], sband["dr"], sband["cut"],
                 sband["frho"], sband["rhor"], sband["z2r"])
    print("wrote s-band file : %s  (preset: %s)" % (sband_path, args.preset))

    dband_name = args.dband
    if args.demo:
        fe_file = args.fe_file
        if fe_file is None:
            here = os.path.dirname(os.path.abspath(__file__))
            fe_file = os.path.join(here, "..", "..", "potentials", "Fe_mm.eam.fs")
        dband = build_demo_dband(fe_file, elements)
        dband_name = args.prefix + ".dband.eam.fs"
        write_eam_fs(dband_name, dband["comments"], dband["elements"],
                     dband["masses"], dband["znums"], dband["a0s"], dband["lats"],
                     dband["nrho"], dband["drho"], dband["nr"], dband["dr"],
                     dband["cut"], dband["frho"], dband["rhor"], dband["z2r"])
        print("wrote demo d-band : %s  (illustrative; both elements are Fe)" % dband_name)
        inp = os.path.join(os.path.dirname(args.prefix),
                           "in." + os.path.basename(args.prefix) + "_twoband")
        write_sample_input(inp, dband_name, sband_path, elements)
        print("wrote sample input: %s" % inp)

    print()
    print("Use in a LAMMPS input with:")
    print("    pair_style hybrid/overlay eam/fs eam/fs")
    if dband_name:
        print("    pair_coeff * * eam/fs 1 %s %s" % (dband_name, " ".join(elements)))
    else:
        print("    pair_coeff * * eam/fs 1 <your-dband>.eam.fs %s" % " ".join(elements))
    print("    pair_coeff * * eam/fs 2 %s %s" % (sband_path, " ".join(elements)))


if __name__ == "__main__":
    main()
