.. index:: velocity/tspin

velocity/tspin command
======================

Syntax
""""""

.. code-block:: LAMMPS

   velocity/tspin group-ID create T seed keyword value ...

* group-ID = ID of group of atoms whose spin velocities will be set
* create = style name of this command
* T = desired spin temperature (temperature units, K in metal units)
* seed = random number seed to use for the velocity generation (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *spinmass* or *mom*

  .. parsed-literal::

       *spinmass* value = ms
         ms = spin mass of each atom, as a multiple of its atomic mass (adim)
       *mom* value = *yes* or *no*
         yes = subtract the mean spin velocity of each atom type
         no = keep the mean spin velocity

Examples
""""""""

.. code-block:: LAMMPS

   velocity/tspin all create 300.0 12345 spinmass 0.0075
   velocity/tspin magnetic create 0.0 12345 spinmass 0.01 mom no

Description
"""""""""""

.. versionadded:: TBD

Assign spin masses and create a random distribution of spin velocities
for inertial spin dynamics.  This is the spin analogue of :doc:`velocity
create <velocity>`; the :doc:`velocity <velocity>` command itself is
unaffected and still sets the atom velocities.

If the *spinmass* keyword is used, the spin mass of every magnetic atom
in the group is set to a multiple of the atomic mass of its type,
:math:`m_s = ms \times m_i`.  Spin masses can equivalently be assigned
with the *spinmass* keyword of :doc:`fix nve/tspin <fix_nve_tspin>`.

The spin velocities of the group are then drawn from a Gaussian
distribution and rescaled so that each atom type separately satisfies

.. math::

   \frac{3}{2} N k_B T = \sum_i \frac{1}{2} m_s \left| \vec{v}^{s}_i \right|^2

If *mom* is *yes*, the mean spin velocity of each type is subtracted
before the rescaling.  Atoms whose spin modulus is zero carry no spin
degrees of freedom and are skipped.

Note: the random number *seed* must be a positive integer.  Each
processor uses the input seed to generate its own unique seed and its
own stream of random numbers, so the generated spin velocities depend
on the number of processors.

Restrictions
""""""""""""

The *velocity/tspin* command is part of the SPIN package.  This command
is only enabled if LAMMPS was built with this package.  See the
:doc:`Build package <Build_package>` page for more info.

This command requires :doc:`atom_style tspin <atom_style>`, and must be
used after the simulation box and the atoms have been defined.

Related commands
""""""""""""""""

:doc:`velocity <velocity>`,
:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`set <set>`

Default
""""""""

The default value for the *mom* keyword is *yes*.  If *spinmass* is not
used, the spin masses are left unchanged.
