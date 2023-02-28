.. index:: fix alchemy

fix alchemy command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID alchemy v_name

* ID, group-ID are documented in :doc:`fix <fix>` command
* alchemy = style name of this fix command
* v_name = variable with name that determines the :math:`\lambda_p` value

Examples
""""""""

.. code-block:: LAMMPS

   fix trans all alchemy v_ramp

Description
"""""""""""

.. versionadded:: TBD

This fix command enables running an "alchemical transformation" MD
simulation between two topologies (i.e. the same number and positions of
atoms, but differences in atom parameters like type, charge, bonds,
angles and so on).  For this a :ref:`multi-partition run <partition>` is
required with exactly two partitions.  During the MD run, the fix will
will determine a factor, :math:`\lambda_p`, for each partition *p* that
will be taken from an equal style or equivalent :doc:`variable
<variable>`.  Typically, this variable would be chose to linearly ramp
*down* from 1.0 to 0.0 for the *first* partition (*p=0*) and linearly
ramp *up* from 0.0 to 1.0 for the *second* partition (*p=1*).  The
forces used for the propagation of the atoms will be the sum of the
forces of the two systems combined and scaled with their respective
:math:`\lambda_p` factor.  This allows to perform transformations that
are not easily possible with :doc:`pair style hybrid/scaled
<pair_hybrid>`, :doc:`fix adapt <fix_adapt>` or :doc:`fix adapt/fep
<fix_adapt_fep>`.

.. note::

   Since the definition of the variable to provide the :math:`\lambda_p` is
   independent in the two partitions, no check is made that the two values
   remain between 0.0 and 1.0 and that they add up to 1.0.  So care needs to
   be taken when defining those variables that this is the case.

Due to the specifics of the implementation, the initial geometry and
dimensions of the system must be exactly the same and the fix will
synchronize them during the run.  It is thus not possible to initialize
the two partitions by reading different data files or creating different
systems from scratch, but rather they have to be started from the same
system and then the desired modifications need to be applied to the
system of the second partition.  The commands :doc:`pair style
hybrid/scaled <pair_hybrid>`, :doc:`fix adapt <fix_adapt>` or :doc:`fix
adapt/fep <fix_adapt_fep>` could be used for simulations where the
requirements for fix alchemy are not given.

The commands below demonstrate how the setup for the second partition
can be done for the example of transforming a pure copper system into a
copper/aluminum bronze.

.. code-block:: LAMMPS

   variable name world pure alloy

   create_box 2 box
   create_atoms 1 box
   pair_style eam/alloy
   pair_coeff * * AlCu.eam.alloy Cu Al

   # replace 5% of copper with aluminum on the second partition only
   variable name world pure alloy
   if "${name} == alloy" then &
     "set type 1 type/fraction 2 0.05 6745234"

   # define ramp variable to combine the two different partitions
   if "${name} == pure" then             &
     "variable ramp equal ramp(1.0,0.0)"    &
   else                                      &
      "variable ramp equal ramp(0.0,1.0)"

   fix 2 all alchemy v_ramp


The ``examples/PACKAGES/alchemy`` folder contains complete example
inputs for this command.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.

This fix stores a global scalar (the current value of :math:`\lambda_p`)
and a global vector of length 3 which contains the potential energy of
the first partition, the second partition and the combined value,
respectively. The global scalar is unitless and "intensive", the vector
is in :doc:`energy units <units>` and "extensive".  These values can be
used by any command that uses a global value from a fix as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the REPLICA package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

There may be only one instance of this fix in use at any time.

This fix requires to perform a :ref:`multi-partition run <partition>`
with *exactly* two partitions.

This fix is *not* compatible with :doc:`load balancing <fix_balance>`.

Related commands
""""""""""""""""

:doc:`compute pressure/alchemy <compute_pressure_alchemy>` command,
:doc:`fix adapt <fix_adapt>` command, :doc:`fix adapt/fep <fix_adapt_fep>`
command, :doc:`pair_style hybrid/scaled <pair_hybrid>` command.

Default
"""""""

none
