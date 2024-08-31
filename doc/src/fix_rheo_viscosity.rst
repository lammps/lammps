.. index:: fix rheo/viscosity

fix rheo/viscosity command
==========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo/viscosity type1 pstyle1 args1 ... typeN pstyleN argsN

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo/viscosity = style name of this fix command
* one or more types and viscosity styles must be appended
* types = lists of types (see below)
* vstyle = *constant* or *power*

  .. parsed-literal::

       *constant* args = *eta*
         *eta* = viscosity

       *power* args = *eta*, *gd0*, *K*, *n*
         *eta* = viscosity
         *gd0* = critical strain rate
         *K* = consistency index
         *n* = power-law exponent

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo/viscosity * constant 1.0
   fix 1 all rheo/viscosity 1 constant 1.0 2 power 0.1 5e-4 0.001 0.5

Description
"""""""""""

.. versionadded:: 29Aug2024

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

The *types* definition is followed by the viscosity style, *vstyle*. Two
options are available, *constant* and *power*. Style *constant* simply
applies a constant value of the viscosity *eta* to each particle of the
assigned type. Style *power* is a Hershchel-Bulkley constitutive equation
for the stress :math:`\tau`

.. math::

   \tau = \left(\frac{\tau_0}{\dot{\gamma}} + K \dot{\gamma}^{n - 1}\right) \dot{\gamma}, \tau \ge \tau_0

where :math:`\dot{\gamma}` is the strain rate and :math:`\tau_0` is the critical
yield stress, below which :math:`\dot{\gamma} = 0.0`. To avoid divergences, this
expression is regularized by defining a critical strain rate *gd0*. If the local
strain rate on a particle falls below this limit, a constant viscosity of *eta*
is assigned. This implies a value of

.. math::
   \tau_0 = \eta \dot{\gamma}_0 - K \dot{\gamma}_0^N


Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with an atom style that includes viscosity
such as atom_style rheo or rheo/thermal. This fix must be used in
conjunction with :doc:`fix rheo <fix_rheo>`. The fix group must be
set to all. Only one instance of fix rheo/viscosity can be defined.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none

