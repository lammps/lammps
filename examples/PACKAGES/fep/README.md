compute fep examples
====================

The examples below illustrate the use of *compute fep*, *fix adapt/fep* and
*pair lj/cut/coul/long/soft* to calculate free energy differences. They present
just one version of many possible ways of constructing a path from the initial
to the final states.

* `CH4hyd` -- free energy of hydration of methane (simple). FEP
  and FDTI methods.

* `SPCEhyd` -- free energy of hydration of SPCE water (simple). FEP
  in separate steps for the LJ sites and the atomic charges.

* `CH4-CF4` -- free energy difference of transforming methane into
  perfluoromethane, in water (quite simple). FEP and BAR methods.

* `CC-CO` -- hydration free energy difference between ethane and methanol
 (a bit more complex). FEP method.
  
* `C7inEthanol` -- insertion / deletion of a hexane molecule in ethanol.
  The hexane is described by the *lj/class2/coul/long/soft* potential.

* `quicktests` -- very short runs with charged Lennard-Jones atoms to test
  *compute fep*, *fix adapt/fep* and *pair lj/cut/coul/long/soft*.

* `ta` -- surface tension of SPCE water without constraints. Test-area method.
