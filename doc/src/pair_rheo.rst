.. index:: pair_style rheo

pair_style rheo command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style rheo cut keyword values

* cut = *quintic* or *CRK0* or *CRK1* or *CRK2*
* zero or more keyword/value pairs may be appended to args
* keyword = *rho/damp* or *artificial/visc*

.. parsed-literal::

     *rho/damp* args = density damping prefactor :math:`\xi` (units?)
     *artificial/visc* args = artificial viscosity prefactor :math:`\zeta` (units?)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style rheo 1.0 quintic rho/damp 1.0 artificial/visc 2.0
   pair_coeff * *

Description
"""""""""""

pair style...

No coefficients are defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples
above.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.  Thus, you need to re-specify the pair_style and
pair_coeff commands in an input script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>` command.  It does not support the *inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`,
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

No density damping or artificial viscous forces are calculated.
