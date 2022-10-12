.. index:: fix numdiff/virial

fix numdiff/virial command
==========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID numdiff/virial Nevery delta

* ID, group-ID are documented in :doc:`fix <fix>` command
* numdiff/virial = style name of this fix command
* Nevery = calculate virial by finite difference every this many timesteps
* delta = magnitude of strain fields (dimensionless)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all numdiff/stress 10 1e-6

Description
"""""""""""

.. versionadded:: 17Feb2022

Calculate the virial stress tensor through a finite difference calculation of
energy versus strain.  These values can be compared to the analytic virial
tensor computed by pair styles, bond styles, etc.  This can be useful for
debugging or other purposes. The specified group must be "all".

This fix applies linear strain fields of magnitude *delta* to all the
atoms relative to a point at the center of the box.  The
strain fields are in six different directions, corresponding to the
six Cartesian components of the stress tensor defined by LAMMPS.
For each direction it applies the strain field in both the positive
and negative senses, and the new energy of the entire system
is calculated after each. The difference in these two energies
divided by two times *delta*, approximates the corresponding
component of the virial stress tensor, after applying
a suitable unit conversion.

.. note::

   It is important to choose a suitable value for delta, the magnitude of
   strains that are used to generate finite difference
   approximations to the exact virial stress.  For typical systems, a value in
   the range of 1 part in 1e5 to 1e6 will be sufficient.
   However, the best value will depend on a multitude of factors
   including the stiffness of the interatomic potential, the thermodynamic
   state of the material being probed, and so on. The only way to be sure
   that you have made a good choice is to do a sensitivity study on a
   representative atomic configuration, sweeping over a wide range of
   values of delta.  If delta is too small, the output values will vary
   erratically due to truncation effects. If delta is increased beyond a
   certain point, the output values will start to vary smoothly with
   delta, due to growing contributions from higher order derivatives. In
   between these two limits, the numerical virial values should be largely
   independent of delta.

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

This fix produces a global vector which can be accessed by various
:doc:`output commands <Howto_output>`, which stores the components of
the virial stress tensor as calculated by finite difference.  The
global vector can only be accessed on timesteps that are multiples
of *Nevery* since that is when the finite difference virial is
calculated. See the examples in *examples/numdiff* directory
to see how this fix can be used to directly compare with
the analytic virial stress tensor computed by LAMMPS.

The order of the virial stress tensor components is *xx*, *yy*, *zz*,
*yz*, *xz*, and *xy*, consistent with Voigt notation. Note that
the vector produced by :doc:`compute pressure <compute_pressure>`
uses a different ordering, with *yz* and *xy* swapped.

The vector values calculated by this compute are
"intensive".  The vector values will be in pressure
:doc:`units <units>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is invoked during :doc:`energy
minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix numdiff <fix_numdiff>`, :doc:`compute pressure <compute_pressure>`

Default
"""""""

none
