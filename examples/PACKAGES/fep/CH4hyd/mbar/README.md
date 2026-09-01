Free Energy of Hydration of Methane with MBAR
=============================================

Example calculation of the free energy of hydration of methane using the
multistate Bennett acceptance ratio (MBAR) method, with LAMMPS *compute mbar*,
*fix adapt/fep* and *pair lj/cut/tip4p/long/soft*.

This is the multistate analogue of the FEP calculation in the parent `CH4hyd`
directory: instead of perturbing one state at a time, *compute mbar* evaluates
the reduced potential of each sampled configuration at every state, producing
one row of the u_kln matrix that MBAR requires.

The transformation is here split into two legs, run in sequence:

* `in-mbar-lj.lmp` -- grow the Lennard-Jones (van der Waals) interactions of
  the methane sites. *fix adapt/fep* holds the soft-core activation parameter
  lambda at 21 equally spaced stages from 0.0 to 1.0 (a window of 50000 steps
  each); *compute mbar* uses the matching grid `0.0 1.0 21`. Reads `data.lmp`,
  writes the reduced potentials to `mbar01-lj.lmp` and the final configuration
  to `data-mbar-lj.lmp`.

* `in-mbar-q.lmp` -- grow the partial charges of the methane sites. *fix
  adapt/fep* holds the charges at 11 stages (window 20000 steps); *compute
  mbar* uses the grids `0.0 -0.24 11` and `0.0 0.06 11`. Reads
  `data-mbar-lj.lmp` (the output of the LJ leg) and writes `mbar01-q.lmp`.

Run the legs in order:

    lmp -in in-mbar-lj.lmp
    lmp -in in-mbar-q.lmp

Each leg writes, every 20 steps, a `fix ave/time ... mode vector` file holding
the instantaneous reduced potentials at every state (no time averaging, so the
raw samples are available for decorrelation). Post-process them with the
scripts in the `tools/fep` directory: `lmp2ukln.py` reshapes the LAMMPS output
into a u_kln array (grouping samples by held state using the window length),
and `mbar.py` runs pymbar (per-state equilibration detection and decorrelation
included) to obtain the free energy difference and a profile along lambda:

    lmp2ukln.py mbar01-lj.lmp 50000 u_kln-lj.npy
    mbar.py 300 u_kln-lj.npy

    lmp2ukln.py mbar01-q.lmp 20000 u_kln-q.npy
    mbar.py 300 u_kln-q.npy

The two contributions (LJ and charge) add up to the free energy of hydration,
dominated by the LJ/cavity term, and can be compared with the FEP result in the
parent directory and with the experimental value of 2.0 kcal/mol.

These example calculations are for tutorial purposes only. The results may not
be of research quality (sampling, lambda spacing, ideal-gas contributions,
etc. are not optimized).
