.. index:: pair_style rheo/react

pair_style rheo/react command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style rheo/react

Examples
""""""""

.. code-block:: LAMMPS

   pair_style rheo/react
   pair_coeff * * 1.0 1.5 1.0 0.05 1.0 100 2.0

Description
"""""""""""

pair style...

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the example above,
or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* :math:`k` (force/distance units)
* :math:`r_max` (distance units)
* :math:`\epsilon` (unitless)
* :math:`\gamma` (force/velocity units)
* :math:`t_form` (time units)
* :math:`r_from_surface` (distance units)

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
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none
