.. index:: compute ptm/atom

compute ptm/atom command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID ptm/atom structures threshold group2-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* ptm/atom = style name of this compute command
* structures = *default* or *all* or any hyphen-separated combination of *fcc*, *hcp*, *bcc*, *ico*, *sc*, *dcub*, *dhex*, or *graphene* = structure types to search for
* threshold = lattice distortion threshold (RMSD)
* group2-ID determines which group is used for neighbor selection (optional, default "all")

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all ptm/atom default 0.1 all
   compute 1 all ptm/atom fcc-hcp-dcub-dhex 0.15 all
   compute 1 all ptm/atom all 0

Description
"""""""""""

Define a computation that determines the local lattice structure
around an atom using the PTM (Polyhedral Template Matching) method.
The PTM method is described in :ref:`(Larsen) <Larsen>`.

Currently, there are seven lattice structures PTM recognizes:

* fcc = 1
* hcp = 2
* bcc = 3
* ico (icosahedral) = 4
* sc (simple cubic) = 5
* dcub (diamond cubic) = 6
* dhex (diamond hexagonal) = 7
* graphene = 8

The value of the PTM structure will be 0 for unknown types and :math:`-1` for
atoms not in the specified compute group.  The choice of structures to search
for can be specified using the "structures" argument, which is a
hyphen-separated list of structure keywords.
Two convenient pre-set options are provided:

* default: fcc-hcp-bcc-ico
* all: fcc-hcp-bcc-ico-sc-dcub-dhex-graphene

The 'default' setting detects the same structures as the Common Neighbor Analysis method.
The 'all' setting searches for all structure types.  A performance penalty is
incurred for the diamond and graphene structures, so it is not recommended to use this option if
it is known that the simulation does not contain these structures.

PTM identifies structures using two steps.  First, a graph isomorphism test is used
to identify potential structure matches.  Next, the deviation is computed between the
local structure (in the simulation) and a template of the ideal lattice structure.
The deviation is calculated as:

.. math::

   \text{RMSD}(\mathbf{u}, \mathbf{v})
    = \min_{s, \mathbf{Q}} \sqrt{\frac{1}{N} \sum\limits_{i=1}^{N}
   {\left\lVert s[\vec{u_i} - \mathbf{\bar{u}}]
                 - \mathbf{Q} \cdot \vec{v_i} \right\rVert}^2}

Here, :math:`\vec u` and :math:`\vec v` contain the coordinates of the local
and ideal structures respectively, :math:`s` is a scale factor, and
:math:`\mathbf Q` is a rotation.  The best match is identified by the lowest
RMSD value, using the optimal scaling, rotation, and correspondence between the
points.

The *threshold* keyword sets an upper limit on the maximum permitted deviation
before a local structure is identified as disordered.  Typical values are in
the range 0.1--0.15, but larger values may be desirable at higher temperatures.
A value of 0 is equivalent to infinity and can be used if no threshold is
desired.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (e.g., each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently or to have multiple compute/dump commands, each with a
*ptm/atom* style. By default the compute processes **all** neighbors
unless the optional *group2-ID* argument is given, then only members
of that group are considered as neighbors.

Output info
"""""""""""

This compute calculates a per-atom array, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

Results are stored in the per-atom array in the following order:

* type
* rmsd
* interatomic distance
* qw
* qx
* qy
* qz

The type is a number from :math:`-1` to 8. The rmsd is a positive real number.
The interatomic distance is computed from the scale factor in the RMSD equation.
The :math:`(qw,qx,qy,qz)` parameters represent the orientation of the local
structure in quaternion form.  The reference coordinates for each template
(from which the orientation is determined) can be found in the
*ptm_constants.h* file in the PTM source directory.
For atoms that are not within the compute group-ID, all values are set to zero.

Restrictions
""""""""""""

This fix is part of the PTM package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute centro/atom <compute_centro_atom>`
:doc:`compute cna/atom <compute_cna_atom>`

Default
"""""""

none

----------

.. _Larsen:

**(Larsen)** Larsen, Schmidt, Schiotz, Modelling Simul Mater Sci Eng, 24, 055007 (2016).
