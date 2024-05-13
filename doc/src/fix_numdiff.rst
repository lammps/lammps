.. index:: fix numdiff

fix numdiff command
===================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID numdiff Nevery delta

* ID, group-ID are documented in :doc:`fix <fix>` command
* numdiff = style name of this fix command
* Nevery = calculate force by finite difference every this many timesteps
* delta = size of atom displacements (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all numdiff 10 1e-6
   fix 1 movegroup numdiff 100 0.01

Description
"""""""""""

Calculate forces through finite difference calculations of energy
versus position.  These forces can be compared to analytic forces
computed by pair styles, bond styles, etc.  This can be useful for
debugging or other purposes.

The group specified with the command means only atoms within the group
have their averages computed.  Results are set to 0.0 for atoms not in
the group.

This fix performs a loop over all atoms in the group.  For each atom
and each component of force it adds *delta* to the position, and
computes the new energy of the entire system.  It then subtracts
*delta* from the original position and again computes the new energy
of the system.  It then restores the original position.  That
component of force is calculated as the difference in energy divided
by two times *delta*.

.. note::

   It is important to choose a suitable value for delta, the magnitude of
   atom displacements that are used to generate finite difference
   approximations to the exact forces.  For typical systems, a value in
   the range of 1 part in 1e4 to 1e5 of the typical separation distance
   between atoms in the liquid or solid state will be sufficient.
   However, the best value will depend on a multitude of factors
   including the stiffness of the interatomic potential, the thermodynamic
   state of the material being probed, and so on. The only way to be sure
   that you have made a good choice is to do a sensitivity study on a
   representative atomic configuration, sweeping over a wide range of
   values of delta.  If delta is too small, the output forces will vary
   erratically due to truncation effects. If delta is increased beyond a
   certain point, the output forces will start to vary smoothly with
   delta, due to growing contributions from higher order derivatives. In
   between these two limits, the numerical force values should be largely
   independent of delta.

.. note::

   The cost of each energy evaluation is essentially the cost of an MD
   timestep.  Thus invoking this fix once for a 3d system has a cost
   of 6N timesteps, where N is the total number of atoms in the system.
   So this fix can be very expensive to use for large systems.
   One expedient alternative is to define the fix for a group containing
   only a few atoms.

----------

The *Nevery* argument specifies on what timesteps the force will
be used calculated by finite difference.

The *delta* argument specifies the size of the displacement each
atom will undergo.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.

This fix produces a per-atom array which can be accessed by various
:doc:`output commands <Howto_output>`, which stores the components of
the force on each atom as calculated by finite difference.  The
per-atom values can only be accessed on timesteps that are multiples
of *Nevery* since that is when the finite difference forces are
calculated. See the examples in *examples/numdiff* directory
to see how this fix can be used to directly compare with
the analytic forces computed by LAMMPS.

The array values calculated by this compute
will be in force :doc:`units <units>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is invoked during :doc:`energy
minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`dynamical_matrix <dynamical_matrix>`, :doc:`fix numdiff/virial <fix_numdiff_virial>`,

Default
"""""""

none
