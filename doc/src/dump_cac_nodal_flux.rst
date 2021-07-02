.. index:: dump cac/nodal/flux

dump cac/nodal/flux command
=============================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/nodal/flux N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* cac/nodal/flux = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = *sizex* *sizey* *sizez* or *sizex* *sizey* *sizez* *dx* *dy* *dz*

Examples
""""""""

.. code-block:: LAMMPS

   dump 1 all cac/flux 1000 CACflux.txt quadrature 2.715 2.715 2.715
   dump 1 all cac/flux 1000 CACflux.txt quadrature 2.715 2.715 2.715 0.2 0.1 0.3

Description
"""""""""""

Periodically outputs the nodal fluxes, computed using the surface integral averaged approach
found in :ref:`(Chen-Diaz) <USER-CAC-Chen-Diaz>`, at the specified interval. In the case of atoms, the algorithm places a user specified
box around the atom and computes the surface fluxes through each plane of the box. In the case of
elements the box is placed around each quadrature point in the case of the *quadrature* option, or
just the quadrature points nearest to element nodes in the case of the *nodal* option. The advantage
of the latter option is that the interpolated flux represents surface defects, or free surfaces,
more accurately whilst neglecting an accurate interpolation of the element interior.
The nodal fluxes are defined as the equivalent nodal fluxes in the case of the *quadrature* option.
The boxes used to compute fluxes are specified using the *sizex*, *sizey*, *sizez*, *dx*, *dy*, and *dz*
arguments to represent the box dimensions and box center displacement from the atom in all three dimension.
If the three offsets are omitted they default to 0, the box is centered around the atom or quadrature point.
Currently, all three offsets must be supplied, even if some are 0, if defining offsets, otherwise an
error is triggered.

The nodal fluxes of all finite elements for each of their respective internal degrees of
freedom will be printed; atom fluxes will also be printed and denoted by the 0 element
type following one line of flux data. Each row will contain 24 values, corresponding to the heat flux
and three stress components on each of the box's bounding planes. The order of surface normals is -x,x,-y,y,-z,z
, with each plane having its four fluxes output in order with heat flux, x stress, y stress, z stress.
Note that the output fluxes are not divided by their respective surface areas,
this is left for the user to decide. A common choice is the specified box surface areas.

----------

NOTE: The :doc:`dump_modify sort <dump_modify>` option
does not work with this dump style.

Restrictions
""""""""""""

This dump style requires a CAC :doc:`atom style <atom_style>`

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump cac/initial/nodes <dump_cac_initial_nodes>`,
:doc:`dump cac/kinetic/energy <dump_cac_kinetic>`, :doc:`dump cac/xyz <dump_cac_xyz>`,
:doc:`dump cac/nodal/positions <dump_cac_nodal_positions>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`

----------

.. _USER-CAC-Chen-Diaz:

**(USER-CAC-Chen-Diaz)** Chen, Youping, and Adrian Diaz. "Physical foundation and consistent formulation of atomic-level fluxes in transport processes." Physical Review E 98.5 (2018): 052113.
