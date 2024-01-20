.. index:: pair_style rheo/solid

pair_style rheo/solid command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style rheo/solid

Examples
""""""""

.. code-block:: LAMMPS

   pair_style rheo/solid
   pair_coeff * * 1.0 1.5 1.0

Description
"""""""""""

pair style...

* :math:`k` (force/distance units)
* :math:`\sigma` (distance units)
* :math:`\gamma` (force/velocity units)

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
:doc:`pair bpm/spring <pair_bpm_spring>`,

Default
"""""""

none
