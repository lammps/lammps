.. index:: pair_style coul/cut/dielectric
.. index:: pair_style coul/long/dielectric
.. index:: pair_style lj/cut/coul/cut/dielectric
.. index:: pair_style lj/cut/coul/cut/dielectric/omp
.. index:: pair_style lj/cut/coul/debye/dielectric
.. index:: pair_style lj/cut/coul/debye/dielectric/omp
.. index:: pair_style lj/cut/coul/long/dielectric
.. index:: pair_style lj/cut/coul/long/dielectric/omp
.. index:: pair_style lj/cut/coul/msm/dielectric
.. index:: pair_style lj/long/coul/long/dielectric

pair_style coul/cut/dielectric command
======================================

pair_style coul/long/dielectric command
=======================================

pair_style lj/cut/coul/cut/dielectric command
=============================================

Accelerator Variants: *lj/cut/coul/cut/dielectric/omp*

pair_style lj/cut/coul/debye/dielectric command
===============================================

Accelerator Variants: *lj/cut/coul/debye/dielectric/omp*

pair_style lj/cut/coul/long/dielectric command
==============================================

Accelerator Variants: *lj/cut/coul/long/dielectric/omp*

pair_style lj/cut/coul/msm/dielectric command
==============================================

pair_style lj/long/coul/long/dielectric command
===============================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/cut/coul/cut/dielectric* or *lj/cut/coul/long/dielectric* or *lj/cut/coul/msm/dielectric* or *lj/long/coul/msm/dielectric*
* args = list of arguments for a particular style

Examples
""""""""

.. code-block:: LAMMPS

   pair_style coul/cut/dielectric 10.0
   pair_coeff * *
   pair_coeff 1 1 9.0

   pair_style lj/cut/coul/cut/dielectric 10.0
   pair_style lj/cut/coul/cut/dielectric 10.0 8.0
   pair_coeff * * 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

   pair_style lj/cut/coul/long/dielectric 10.0
   pair_style lj/cut/coul/long/dielectric 10.0 8.0
   pair_coeff * * 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

Used in input scripts:

   .. parsed-literal::

      examples/PACKAGES/dielectric/in.confined
      examples/PACKAGES/dielectric/in.nopbc

Description
"""""""""""

All these pair styles are derived from the corresponding pair styles
without the *dielectric*\ suffix. In addition to computing atom forces
and energies, these pair styles compute the electrical field vector
at each atom, which are to be used in the :doc:`fix polarize <fix_polarize>` commands.

These pair styles should be used with :doc:`atom_style dielectric <atom_style>`,
which uses atom charges rescaled by their local dielectric constant.

The styles lj/cut/coul/long/dielectric, lj/cut/coul/msm/dielectric, and
lj/long/coul/long/dielectric should be used with their kspace style counterparts,
namely, pppm/dielectric, pppm/disp/dielectric, and msm/dielectric, respectively.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distances for this pair style can be mixed.  The default
mix value is *geometric*\ .  See the "pair_modify" command for details.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

These styles are part of the DIELECTRIC package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix polarize <fix_polarize>`, :doc:`read_data <read_data>`

Default
"""""""

none

