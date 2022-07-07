.. index:: compute pace/extrapolation

compute pace/extrapolation command
==================================

Syntax
""""""

.. parsed-literal::

   compute ID all pace/extrapolation

* ID is documented in :doc:`compute <compute>` command
* pace/extrapolation = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute pace_gamma all pace/extrapolation

Description
"""""""""""

Define a computation that calculates both per-atom and per-structure extrapolation grades for PACE interatomic potential.
Pair style :doc:`pair_style pace/extrapolation <pair_pace>` must be instantiated before.
Extrapolation grades are computed by `pair_style pace/extrapolation` every *gamma_freq* steps (see :doc:`pair_style pace/extrapolation <pair_pace>`),
but `compute pace/extrapolation` will invoke extra calculations with this pair style if necessary.
For better performance, it is recommended to use the same values of *gamma_freq* and
the frequency of compute style callers, i.e. `dump` or `thermo`.

.. code-block:: LAMMPS

   compute pace_gamma all pace/extrapolation

   # per-structure extrapolation grade c_pace_gamma
   thermo_style custom step etotal temp press c_pace_gamma

   # per-atom extrapolation grade c_pace_gamma
   dump 1 all custom 100 my.dump id type mass x y z c_pace_gamma


Output info
"""""""""""

This compute calculates both per-atom vector and per-structure scalar,
which can be accessed by any command that uses per-atom and/or per-structure values from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

All values are unitless.

Restrictions
""""""""""""

 Pair style :doc:`pair_style pace/extrapolation <pair_pace>` must be instantiated before.

 group-ID always corresponds to the group atoms used by `pair_style pace/extrapolation` and by default is `all`.

Related commands
""""""""""""""""

:doc:`pair_style pace/extrapolation <pair_pace>`, :doc:`dump custom <dump>`, :doc:`thermo custom <thermo>`

Default
"""""""

`compute pace_gamma all pace/extrapolation`
