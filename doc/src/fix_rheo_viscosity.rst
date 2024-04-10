.. index:: fix rheo/viscosity

fix rheo/viscosity command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo/viscosity type1 pstyle1 args1 ... typeN pstyleN argsN

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo/viscosity = style name of this fix command
* one or more types and viscosity styles must be appended
* types = lists of types (see below)
* vstyle = *constant*

  .. parsed-literal::

       *constant* args = *eta*
         *eta* = viscosity

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo/viscosity * constant 1.0
   fix 1 all rheo/viscosity 1 constant 1.0 2 constant 2.0

Description
"""""""""""

.. versionadded:: TBD

This fix defines a viscosity for RHEO particles. One can define different
viscosities for different atom types, but a viscosity must be specified for
every atom type.

One first defines the atom *types*. A wild-card asterisk can be used in place
of or in conjunction with the *types* argument to set the coefficients for
multiple pairs of atom types.  This takes the form "\*" or "\*n" or "m\*"
or "m\*n".  If :math:`N` is the number of atom types, then an asterisk with
no numeric values means all types from 1 to :math:`N`.  A leading asterisk
means all types from 1 to n (inclusive).  A trailing asterisk means all types
from m to :math:`N` (inclusive).  A middle asterisk means all types from m to n
(inclusive).

The *types* definition is followed by the viscosity style, *vstyle*. Currently,
the only option is *constant*. Style *constant* simply applies a constant value
of the viscosity *eta* to each particle of the assigned type.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with an atom style that includes viscosity
such as atom_style rheo or rheo/thermal. This fix must be used in
conjuction with :doc:`fix rheo <fix_rheo>`. The fix group must be
set to all. Only one instance of fix rheo/viscosity can be defined.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none
