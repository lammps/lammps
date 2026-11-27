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
raw samples are available for decorrelation), with 15 significant digits
(`format " %.15g"`) since the default format of 6 digits loses accuracy.
Post-process them with the scripts in the `tools/fep` directory:
`lmp2ukln.py` reshapes the LAMMPS output
into a u_kln array (grouping samples by held state using the window length),
and `mbar.py` runs pymbar (per-state equilibration detection and decorrelation
included) to obtain the free energy difference, here in kcal/mol (`real`
units), and a profile along the states:

    lmp2ukln.py mbar01-lj.lmp 50000 u_kln-lj.npy
    mbar.py real 300 u_kln-lj.npy -g 0.0 1.0

    lmp2ukln.py mbar01-q.lmp 20000 u_kln-q.npy
    mbar.py real 300 u_kln-q.npy

The log files of the two legs, `log.mbar-lj` and `log.mbar-q` (run on 8 MPI
processes), are provided for comparison. The output files with the reduced
potentials are too large to be included; post-processing them as above gave:

    LJ leg:      3.93 +/- 0.22 kT  =  2.34  +/- 0.13  kcal/mol
    charge leg:  0.02 +/- 0.01 kT  =  0.010 +/- 0.007 kcal/mol

The numbers of a different run will differ from these within the statistical
uncertainty, since the trajectories depend on the number of processes and on
the platform.

The two contributions (LJ and charge) add up to the free energy of hydration,
2.35 kcal/mol, dominated by the LJ/cavity term. This can be compared with the
FEP result in the parent directory (2.12 kcal/mol from `fep01`), with the
literature value for these force field models, 2.27 kcal/mol, and with the
experimental value of 2.0 kcal/mol.

These example calculations are for tutorial purposes only. The results may not
be of research quality (sampling, lambda spacing, ideal-gas contributions,
etc. are not optimized).
