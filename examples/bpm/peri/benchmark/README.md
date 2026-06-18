# bpm/peri vs PERI performance benchmark

These decks compare the per-step force-computation cost of the BPM peridynamics
styles (`pair_style` + `bond_style bpm/peri`) against the legacy PERI package
(`pair_style peri/*`) on an identical system.

## Setup

- Periodic 16x16x16 simple-cubic block: 4096 nodes, lattice spacing 0.0005,
  horizon 0.0015001 (3 lattice spacings) -> 122 bonds per interior node.
- Run from the reference configuration (zero net force) for 200 steps, so the
  measured cost is the steady per-step force loop without dynamics, fracture,
  or neighbor-list rebuilds.
- Material: K = 14.9 GPa, density 2200, nodal volume 1.25e-10 (PMB uses the
  derived micromodulus c = 1.6863e22); break threshold s00 = 0.5 (never reached
  at the reference, so no bonds break).
- `in.bench.legacy.<model>` uses `atom_style peri`; `in.bench.bpm.<model>` uses
  `atom_style bond` + `fix property/atom d_vfrac` + `bond_style bpm/peri`.

To benchmark a different model, edit the `pair_style`/`pair_coeff` (legacy) or
`bond_coeff` (BPM) line:

| model | legacy `pair_coeff * *`                                  | BPM `bond_coeff 1`                                |
|-------|----------------------------------------------------------|---------------------------------------------------|
| pmb   | `1.6863e22 0.0015001 0.5 0.25`                           | `pmb 1.6863e22 0.0015001 0.5 0.25`                |
| lps   | `14.9e9 14.9e9 0.0015001 0.5 0.25`                      | `lps 14.9e9 14.9e9 0.0015001 0.5 0.25`            |
| ves   | `14.9e9 14.9e9 0.0015001 0.5 0.25 0.5 0.001`            | `ves 14.9e9 14.9e9 0.5 0.001 0.0015001 0.5 0.25`  |
| eps   | `14.9e9 14.9e9 0.0015001 0.5 0.25 10.0e8`               | `eps 14.9e9 14.9e9 10.0e8 0.0015001 0.5 0.25`     |

## Results (single MPI rank, 200 steps, "Loop time" in seconds)

| model | legacy peri/* | bpm/peri | speedup |
|-------|---------------|----------|---------|
| pmb   | 2.95          | 0.82     | 3.6x    |
| lps   | 7.00          | 1.30     | 5.4x    |
| ves   | 7.96          | 1.80     | 4.4x    |
| eps   | 10.34         | 1.74     | 5.9x    |

np=4 (LPS): legacy 1.92 s, bpm/peri 0.46 s -> 4.1x.

## Why BPM is faster

The legacy PERI pair styles walk a per-node partner list that is **double
counted** (each bond is visited from both endpoints, ~244 visits/node) and
compute the bond force and the short-range contact in the same pass.  The BPM
implementation walks the bond list **once** (single visit, ~122/node, tag
ordered) and the contact `pair_style bpm/peri` has its own short-range neighbor
list, which on a fully bonded lattice is almost empty (all near pairs are
excluded as bonded).  So BPM does roughly half the bond work plus a nearly free
contact pass.  The breakdown (BPM, pmb): bond 0.69 s, pair 0.11 s, comm 0.02 s
-- the forward communication of the dilatation/break state for the state-based
models is a negligible fraction, not the bottleneck.

Hardware/build will shift the absolute numbers; the relative ordering
(bpm/peri several times faster than legacy peri/*) is the robust observation.
