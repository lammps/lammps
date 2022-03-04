Free Energy of Hydration of Methane
===================================

Example calculation of the free energy of hydration of methane with
LAMMPS using *compute fep*, *fix adapt/fep* and *pair lj/cut/coul/long/soft*.

Methane is represented by the 5-site OPLS-AA model (1 molecule). Water
is represented by the 4-site TIP4P-Ew model (360 molecules). Interactions
of sites that are being created or deleted are treated using Lennard-Jones
and Coulomb potentials  with soft cores (*pair lj/cut/coul/long/soft*)
in order to avoid singularities. Trajectories are at constant  NpT, so
corrections for the fluctuating volume are included. 

The following directories contain input files and results for
calculations using free-energy perturbation (FEP) and
finite-difference thermodynamic integration (FDTI):

* `mols` -- molecule description files and force field database used
  to create the initial configuration used for the simulations
  `data.lmp`

* `fep01` -- Calculation using FEP, multi-stage creation of a methane
  molecule. Results in `fep01.fep`

* `fep10` -- Calculation using FEP, multi-stage deletion of a methane
  molecule. Results in `fep10.fep`

* `fdti01` -- Calculation using TI/FDTI, creation of a methane
  molecule. Results in `fdti01.fep`

* `fdti10` -- Calculation using TI/FDTI, deletion a methane
  molecule. Results in `fdti10.fep`

The Python scripts `fep.py`, `fdti.py` and `fti.py` found in the
`tools` directory can be used to calculate the free-energy differences
corresponding to the transformations above:

    fep.py 300 < fep01.fep

    fep.py 300 < fep10.fep

    fdti.py 300 0.002 < fdti01.fep

    fdti.py 300 0.002 < fdti10.fep

    nti.py 300 0.002 < fdti01.fep

    nti.py 300 0.002 < fdti10.fep

The outputs are in kcal/mol and can be compared with the experimental
value of 2.0 kcal/mol, or with a simulation value from the literature
obtained with the same force field models used here: 2.27 kcal/mol
[MR Shirts, VS Pande, J Chem Phys 122 (2005) 134508](https://doi.org/10.1063/1.1877132).

These example calculations are for tutorial purposes only. The results
may not be of research quality (not enough sampling, size of the step
in lambda or of the delta for numerical derivative not optimized, no
evaluation of ideal-gas contributions, etc.)
