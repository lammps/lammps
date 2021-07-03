Free Energy of Hydration of SPCE Water
======================================

Example calculation of the free energy of hydration of water with
LAMMPS using *compute fep*, *fix adapt/fep* and *pair lj/cut/soft*.

The Lennard-Jones sites and the electrostatic charges are
created/annihilated in separate runs, which simplifies the use of
*fix adapt/fep* and *compute fep*. The Lennard-Jones sites are handled
using soft core potentials (*pair lj/cut/soft*). Trajectories are at
constant  NpT, so corrections for the fluctuating volume are included. 

The following directories contain input files and results for
calculations using free-energy perturbation (FEP):

* `mols` -- molecule description files and force field database used
  to create the initial configuration used for the simulations
  `data.lmp`

* `fep01` -- Calculation using FEP, multi-stage creation of one SPC/E
  molecule, LJ and q. Results in `fep01-lj.fep` and `fep01-lj.fep`

* `fep10` -- Calculation using FEP, multi-stage deletion of one SPC/E
  molecule, q and LJ. Results in `fep10-q.fep` and `fep10-lj.fep`

The Python script `fep.py` found in the
`tools` directory can be used to calculate the free-energy differences
corresponding to the transformations above:

    fep.py 300 < fep01-lj.fep

    fep.py 300 < fep01-q.fep

    fep.py 300 < fep10-q.fep

    fep.py 300 < fep10-lj.fep

The outputs are in kcal/mol and can be compared with the experimental
value of -6.3 kcal/mol, or with a simulation value from the literature
of -6.7 kcal/mol
[Gonçalves, Stassen, Pure Appl Chem 76 (2004) 231](https://doi.org/10.1351/pac200476010231).

These example calculations are for tutorial purposes only. The results
may not be of research quality (not enough sampling, size of the step
in lambda not optimized, etc.)
