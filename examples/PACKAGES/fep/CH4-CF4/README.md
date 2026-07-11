Methane to Tetrafluoromethane in Water
======================================

Example calculation of the difference in free energy of hydration upon
transforming methane into tetrafluoromethane with LAMMPS using
*compute fep*, *fix adapt/fep* and *pair lj/cut/coul/long/soft*.

Methane and tetrafluoromethane are represented by the OPLS-AA force
field (1 molecule). Water is represented by the 3-site SPC/E model
(360 molecules).

The procedure used to perform the alchemical transformation is the
following:

* The dual topology approach is used, therefore all the atoms of
  methane and perfluorommethane are present throughout the simulation,
  only some of them are dummy sites at the endpoints of the
  transformation. Masses and intramolecular terms (bonds, angles, dihedrals)
  are not changed.

* Interactions of sites that are being created (from dummy sites, XD-) or
  deleted (to become dummy sites X-D) are treated using soft-core versions
  of the Lennard-Jones and Coulomb potentials (*pair
  lj/cut/coul/long/soft*) in order to avoid singularities. The
  exponent of the coupling parameter lambda in the soft-core pair
  potentials was in this example n = 1.

* Atoms of the fragments being created should not interact with atoms being
  deleted. In this small solute this is guaranteed by the exclusion of 1-2
  and 1-3 interactions, so no modifications of the *pair_coeff* are needed.
  Otherwise the LJ epsilon of those interactions should be zeroed.

* Eletrostatic charges that are modified are varied linearly from the
  initial to the final values. This keeps the overall charge of the
  molecule constant, which is good for the long range electrostatics
  (the coupling parameter lambda has no effect on the kspace terms).

The following directories contain input files and results for
calculations using the free-energy perturbation and Bennet's acceptance
ratio (BAR) method:

* `mols` -- molecule description files and force field database used
  to create the initial configurations used for the simulations
  `data.0.lmp` and `data.1.lmp`

* `fep01` -- Calculation using FEP, 20-step transformation of a CH4
  molecule into CF4, constant NpT. Results in `fep01.fep`

* `fep10` -- Calculation using FEP, 20-step transformation of a
  CF4 molecule into CH4, constant NpT. Results in `fep10.fep`

* `bar01` -- Calculation using BAR, 1-step transformation of a CH4
  molecule into CF4, constant NVT, Results in `bar01.fep`

* `bar10` -- Calculation using BAR, 1-step transformation of a
  CF4 molecule into CH4, constant NVT. Results in `bar10.fep`

The Python scripts `fep.py` and  `bar.py` found in the `tools` directory
can be used to calculate the free-energy difference corresponding to the
transformation:

    fep.py 300 < fep01.fep

    fep.py 300 < fep10.fep

    bar.py 300 bar01.lmp bar10.lmp

The outputs are in kcal/mol and can be compared with the experimental
value of 1.2 kcal/mol, and also with a simulation value from the literature
using a different force field): 0.8 kcal/mol
[Gough, Pearlman, Kolmann, J Chem Phys 99 (1993) 9103](https://doi.org/10.1063/1.465525).
This is a small free energy difference so consider the absolute discrepancies.

These example calculations are for tutorial purposes only. The results
may not be of research quality (not enough sampling, size of the step
in lambda or of the delta for numerical derivative not optimized, no
evaluation of ideal-gas contributions, etc.)
