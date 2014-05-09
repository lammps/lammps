Methane to Tetrafluoromethane in Water
======================================

Example calculation of the difference in free energy of hydration upon
transforming methane into tetrafluoromethane with LAMMPS using
*compute fep* and *fix adapt*.

Methane and tetrafluoromethane are represented by the OPLS-AA force
field (1 molecule). Water is represented by the 3-site SPC/E model
(360 molecules).

The strategy used to perform the alchemical transformation is the
following:

* The dual topology approach is used, therefore all the atoms of
  methane and perfluorommethane are present throughout the simulation,
  only some of them are dummy sites at the endpoints of the
  transformation. Masses and intramolecular terms (bond lengths,
  angles, dihedrals) are not changed.

* Interactions of sites that are being created (from dummy sites) or
  deleted (to become dummy sites) are treated using soft-core verions
  of the Lennard-Jones and Coulomb potentials (*pair
  lj/cut/coul/long/soft*) in order to avoid singularities. The
  exponent of the coupling parameter lambda in the soft-core pair
  potentials was in this example n = 1.

* Eletrostatic charges that are modified are varied linearly from the
  initial to the final values. This keeps the overall charge of the
  molecule constant, which is good for the long range electrostatics
  (the coupling parameter lambda has no effect on the kspace terms).

The following directories contain input files and results for
calculations using Bennet's acceptance ratio (BAR) method:

* `bar01` -- Calculation using BAR, 1-step transformation of a CH4
  molecule into CF4. Results in `bar01.lmp`

* `bar10` -- Calculation using BAR, 1-step transformation of a
  CF4 molecule into CH4. Results in `bar10.lmp`

The Python script `bar.py` found in the `tools` directory can
be used to calculate the free-energy difference corresponding to the
transformation:

    bar.py 300 bar01.lmp bar10.lmp

The outputs are in kcal/mol and can be compared with the experimental
value of 1.2 kcal/mol and with a simulation value from the literature
(using a different force field): 0.8 kcal/mol
[Gough, Pearlman, Kolmann, J Chem Phys 99 (1993) 9103]. This small
free energy difference is difficult to predict.

These example calculations are for tutorial purposes only. The results
may not be of research quality (not enough sampling, size of the step
in lambda or of the delta for numerical derivative not optimized, no
evaluation of ideal-gas contributions, etc.)
