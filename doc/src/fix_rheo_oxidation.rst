.. index:: fix rheo/oxidation

fix rheo/oxidation command
==========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo/oxidation cut btype

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo/oxidation = style name of this fix command
* cut = maximum bond length (distance units)
* btype = type of bonds created

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo/oxidation 1.5 2

Description
"""""""""""

This fix...

Each list consists of a series of type
ranges separated by commas. The range can be specified as a
single numeric value, or a wildcard asterisk can be used to specify a range
of values.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  For
example, if M = the number of atom types, then an asterisk with no numeric
values means all types from 1 to M.  A leading asterisk means all types
from 1 to n (inclusive).  A trailing asterisk means all types from n to M
(inclusive).  A middle asterisk means all types from m to n (inclusive).
Note that all atom types must be included in exactly one of the N collections.

While the *Tfreeze* keyword is optional, the *conductivity* and
*specific/heat* keywords are mandatory.

Multiple instances of this fix may be defined to apply different
properties to different groups. However, the union of fix groups
across all instances of fix rheo/thermal must cover all atoms.
If there are multiple instances of this fix, any intersections in
the fix groups will lead to incorrect thermal integration.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with an atom style that includes temperature,
heatflow, and conductivity such as atom_tyle rheo/thermal This fix
must be used in conjuction with :doc:`fix rheo <fix_rheo>` with the
*thermal* setting.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none
