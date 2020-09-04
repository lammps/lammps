.. index:: dump cac/atom/flux

dump cac/atom/flux command
=====================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/atom/flux N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms/elements to be dumped
* cac/atom/flux = style of dump command
* N = dump every this many timesteps
* file = name of file to write dump info to

Examples
""""""""

.. code-block:: LAMMPS

   dump mesh all cac/atom/flux 1000 CACflux.txt
   dump mesh core cac/atom/flux 50 CACmesh.txt

Description
"""""""""""

Periodically outputs the positions and surface integral averaged fluxes of a set of
virtual atoms and real atoms for the specified group. See :doc:`dump_cac_nodal_flux <dump_cac_nodal_flux>`
for a description of how the surface fluxes are computed; note the output is
not surface normalized. These virtual atoms are the result of converting a
CAC element into atomistic resolution through the interpolation fields. The purpose of
this is to print small regions of a CAC model that may require the use of
an atomistic diagnostic tool such as OVITO. The tag associated with virtual atoms is that
of the CAC element. The *charge* keyword can be used to output per atom charges for atom
styles that define charge.

.. warning::

   This dump requires :doc:`dump_cac_nodal_flux <dump_cac_nodal_flux>` invoked before it.

.. warning::

   CAC elements can be quite large, as a result a conversion to atomistic
   resolution may result in a massive amount of data if your group size is not carefully
   chosen. Currently you cannot elect to output a piece of a CAC element as virtual atoms,
   if an element is included in the group (because its centroid is in the group region etc.)
   then all of its material will be converted to virtual atoms and output.

----------

NOTE: The :doc:`dump_modify sort <dump_modify>` option
does not work with this dump style.
Restrictions
""""""""""""

This dump style requires a CAC :doc:`atom style <atom_style>`

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump cac/initial/nodes <dump_cac_initial_nodes>`,
:doc:`dump cac/kinetic/energy <dump_cac_kinetic>`, :doc:`dump cac/nodal/positions <dump_cac_nodal_positions>`,
:doc:`dump cac/nodal/velocities <dump_cac_nodal_velocities>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
