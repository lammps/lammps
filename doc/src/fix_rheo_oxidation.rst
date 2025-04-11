.. index:: fix rheo/oxidation

fix rheo/oxidation command
==========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo/oxidation cut btype rsurf

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo/oxidation = style name of this fix command
* cut = maximum bond length (distance units)
* btype = type of bonds created
* rsurf = distance from surface to create bonds (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo/oxidation 1.5 2 0.0
   fix 1 all rheo/oxidation 1.0 1 2.0

Description
"""""""""""

.. versionadded:: 29Aug2024

This fix dynamically creates bonds on the surface of fluids to
represent physical processes such as oxidation. It is intended
for use with bond style :doc:`bond rheo/shell <bond_rheo_shell>`.

Every timestep, particles check neighbors within a distance of *cut*.
This distance must be smaller than the kernel length defined in
:doc:`fix rheo <fix_rheo>`. Bonds of type *btype* are created between
a fluid particle and either a fluid or solid neighbor. The fluid particles
must also be on the fluid surface, or within a distance of *rsurf* from
the surface. This process is further described in
:ref:`(Clemmer) <howto_rheo_clemmer2>`.

If used in conjunction with solid bodies, such as those generated
by the *react* option of :doc:`fix rheo/thermal <fix_rheo_thermal>`,
it is recommended to use a :doc:`hybrid bond style <bond_hybrid>`
with different bond types for solid and oxide bonds.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with the bond style :doc:`rheo/shell <bond_rheo_shell>`
and :doc:`fix rheo <fix_rheo>` with surface detection enabled.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`bond rheo/shell <bond_rheo_shell>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

none

----------

.. _howto_rheo_clemmer2:

**(Clemmer)** Clemmer, Pierce, O'Connor, Nevins, Jones, Lechman, Tencer, Appl. Math. Model., 130, 310-326 (2024).
