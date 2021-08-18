Ethane to Methanol in Water
===========================

Example calculation of the difference in free energy of hydration upon
transforming ethane into methanol with LAMMPS using *compute fep*,
*fix adapt/fep* and *pair lj/cut/coul/long/soft*.

Ethane and methanol are represented by the OPLS-AA force field (1
molecule). Water is represented by the 3-site SPC/E model (360
molecules).

The procedure used to perform the alchemical transformation is the
following:

* The dual topology approach is used, therefore all the atoms of
  ethane and methanol are present throughout the simulation, only some
  of them are dummy sites at the endpoints of the
  transformation. Masses and intramolecular terms (bond lengths,
  angles, dihedrals) are not changed.

* Interactions of sites that are being created (from dummy sites, XD-) or
  deleted (to become dummy sites, X-D) are treated using soft-core verions
  of the Lennard-Jones and Coulomb potentials (*pair
  lj/cut/coul/long/soft*) in order to avoid singularities. The
  exponent of the coupling parameter lambda in the soft-core pair
  potentials was in this example n = 1.

* In order to avoid catastrophic overlaps between the fragments being created
  and being deleted, interactions between certain atoms need to be zeroed.
  In this case this concerns the H atoms of the methyl group and the O atom
  of the OH group involved in the transformation. Other pairs of atoms
  belonging to these fragments are either within the 1-2 and 1-3 neighbor
  relations, or are H atoms with no LJ site (H of the OH group).

* Eletrostatic charges that are modified are varied linearly from the
  initial to the final values. This keeps the overall charge of the
  molecule constant, which is good for the long range electrostatics
  (the coupling parameter lambda has no effect on the kspace terms).

The following directories contain input files and results for
calculations using free-energy perturbation (FEP), thermodynamic
integration (TI/FDTI) and Bennet's acceptance ratio methods:

* `mols` -- Molcule description files and force field database used to
  create the initial configurations for the simulations `data.0.lmp`
  and `data.1.lmp`

* `fep01` -- Calculation using FEP, multi-stage transformation of an
  ethane molecule into methanol. Results in `fep01.lmp`

* `fep10` -- Calculation using FEP, multi-stage transformation of a
  methanol molecule into ethane. Results in `fep10.lmp`

The free-energy profiles can be observed by plotting the values in the
third column of the results files. The Python script `fep.py`, found in
the `tools` directory, can be used to calculate the free-energy differences
corresponding to the above transformations:

    fep.py 300 < fep01.lmp

    fep.py 300 < fep10.lmp

The outputs are in kcal/mol and can be compared with the experimental
value of -6.93 kcal/mol and with simulation value from the literature
(obtained with different force field parameters):
-6.7 kcal/mol [Jorgensen, Ravimohan, J Chem Phys 83 (1985) 3050](https://doi.org/10.1063/1.449208),
-6.8 kcal/mol [Goette, Grubmüller, J Comp Chem 30 (2007) 447](https://doi.org/10.1002/jcc.21073).

These example calculations are for tutorial purposes only. The results
may not be of research quality (not enough sampling, size of the step
in lambda or of the delta for numerical derivative not optimized, no
evaluation of ideal-gas contributions, etc.)
