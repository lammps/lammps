.. index:: compute voronoi/atom

compute voronoi/atom command
============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID voronoi/atom keyword arg ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* voronoi/atom = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *only_group* or *surface* or *radius* or *edge_histo* or *edge_threshold*
  or *face_threshold* or *neighbors* or *peratom*

  .. parsed-literal::

       *only_group* = no arg
       *occupation* = no arg
       *surface* arg = sgroup-ID
         sgroup-ID = compute the dividing surface between group-ID and sgroup-ID
           this keyword adds a third column to the compute output
       *radius* arg = v_r
         v_r = radius atom style variable for a poly-disperse Voronoi tessellation
       *edge_histo* arg = maxedge
         maxedge = maximum number of Voronoi cell edges to be accounted in the histogram
       *edge_threshold* arg = minlength
         minlength = minimum length for an edge to be counted
       *face_threshold* arg = minarea
         minarea = minimum area for a face to be counted
       *neighbors* value = *yes* or *no* = store list of all neighbors or no
       *peratom* value = *yes* or *no* = per-atom quantities accessible or no

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all voronoi/atom
   compute 2 precipitate voronoi/atom surface matrix
   compute 3b precipitate voronoi/atom radius v_r
   compute 4 solute voronoi/atom only_group
   compute 5 defects voronoi/atom occupation
   compute 6 all voronoi/atom neighbors yes

Description
"""""""""""

Define a computation that calculates the Voronoi tessellation of the
atoms in the simulation box.  The tessellation is calculated using all
atoms in the simulation, but non-zero values are only stored for atoms
in the group.

By default two per-atom quantities are calculated by this compute.
The first is the volume of the Voronoi cell around each atom.  Any
point in an atom's Voronoi cell is closer to that atom than any other.
The second is the number of faces of the Voronoi cell. This is
equal to the number of nearest neighbors of the central atom,
plus any exterior faces (see note below). If the *peratom* keyword
is set to "no", the per-atom quantities are still calculated,
but they are not accessible.

----------

If the *only_group* keyword is specified the tessellation is performed
only with respect to the atoms contained in the compute group. This is
equivalent to deleting all atoms not contained in the group prior to
evaluating the tessellation.

If the *surface* keyword is specified a third quantity per atom is
computed: the Voronoi cell surface of the given atom. *surface* takes
a group ID as an argument. If a group other than *all* is specified,
only the Voronoi cell facets facing a neighbor atom from the specified
group are counted towards the surface area.

In the example above, a precipitate embedded in a matrix, only atoms
at the surface of the precipitate will have non-zero surface area, and
only the outward facing facets of the Voronoi cells are counted (the
hull of the precipitate). The total surface area of the precipitate
can be obtained by running a "reduce sum" compute on c_2[3]

If the *radius* keyword is specified with an atom style variable as
the argument, a poly-disperse Voronoi tessellation is
performed. Examples for radius variables are

.. code-block:: LAMMPS

   variable r1 atom (type==1)*0.1+(type==2)*0.4
   compute radius all property/atom radius
   variable r2 atom c_radius

Here v_r1 specifies a per-type radius of 0.1 units for type 1 atoms
and 0.4 units for type 2 atoms, and v_r2 accesses the radius property
present in atom_style sphere for granular models.

The *edge_histo* keyword activates the compilation of a histogram of
number of edges on the faces of the Voronoi cells in the compute
group. The argument *maxedge* of the this keyword is the largest number
of edges on a single Voronoi cell face expected to occur in the
sample. This keyword adds the generation of a global vector with
*maxedge*\ +1 entries. The last entry in the vector contains the number of
faces with more than *maxedge* edges. Since the polygon with the
smallest amount of edges is a triangle, entries 1 and 2 of the vector
will always be zero.

The *edge_threshold* and *face_threshold* keywords allow the
suppression of edges below a given minimum length and faces below a
given minimum area. Ultra short edges and ultra small faces can occur
as artifacts of the Voronoi tessellation. These keywords will affect
the neighbor count and edge histogram outputs.

If the *occupation* keyword is specified the tessellation is only
performed for the first invocation of the compute and then stored.
For all following invocations of the compute the number of atoms in
each Voronoi cell in the stored tessellation is counted. In this mode
the compute returns a per-atom array with 2 columns. The first column
is the number of atoms currently in the Voronoi volume defined by this
atom at the time of the first invocation of the compute (note that the
atom may have moved significantly). The second column contains the
total number of atoms sharing the Voronoi cell of the stored
tessellation at the location of the current atom. Numbers in column
one can be any positive integer including zero, while column two
values will always be greater than zero. Column one data can be used
to locate vacancies (the coordinates are given by the atom coordinates
at the time step when the compute was first invoked), while column two
data can be used to identify interstitial atoms.

If the *neighbors* value is set to yes, then this compute creates a
local array with 3 columns. There is one row for each face of each
Voronoi cell. The 3 columns are the atom ID of the atom that owns the
cell, the atom ID of the atom in the neighboring cell (or zero if the
face is external), and the area of the face.  The array can be
accessed by any command that uses local values from a compute as
input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options. More specifically, the array can be
accessed by a :doc:`dump local <dump>` command to write a file
containing all the Voronoi neighbors in a system:

.. code-block:: LAMMPS

   compute 6 all voronoi/atom neighbors yes
   dump d2 all local 1 dump.neighbors index c_6[1] c_6[2] c_6[3]

If the *face_threshold* keyword is used, then only faces
with areas greater than the threshold are stored.

----------

The Voronoi calculation is performed by the freely available `Voro++ package <voronoi_>`_, written by Chris Rycroft at UC Berkeley and LBL,
which must be installed on your system when building LAMMPS for use
with this compute.  See instructions on obtaining and installing the
Voro++ software in the src/VORONOI/README file.

.. _voronoi: http://math.lbl.gov/voro++/

.. note::

   The calculation of Voronoi volumes is performed by each
   processor for the atoms it owns, and includes the effect of ghost
   atoms stored by the processor.  This assumes that the Voronoi cells of
   owned atoms are not affected by atoms beyond the ghost atom cut-off
   distance.  This is usually a good assumption for liquid and solid
   systems, but may lead to underestimation of Voronoi volumes in low
   density systems.  By default, the set of ghost atoms stored by each
   processor is determined by the cutoff used for
   :doc:`pair_style <pair_style>` interactions.  The cutoff can be set
   explicitly via the :doc:`comm_modify cutoff <comm_modify>` command.  The
   Voronoi cells for atoms adjacent to empty regions will extend into
   those regions up to the communication cutoff in x, y, or z.  In that
   situation, an exterior face is created at the cutoff distance normal
   to the x, y, or z direction.  For triclinic systems, the exterior face
   is parallel to the corresponding reciprocal lattice vector.

.. note::

   The Voro++ package performs its calculation in 3d.  This will
   still work for a 2d LAMMPS simulation, provided all the atoms have the
   same z coordinate. The Voronoi cell of each atom will be a columnar
   polyhedron with constant cross-sectional area along the z direction
   and two exterior faces at the top and bottom of the simulation box. If
   the atoms do not all have the same z coordinate, then the columnar
   cells will be accordingly distorted. The cross-sectional area of each
   Voronoi cell can be obtained by dividing its volume by the z extent of
   the simulation box.  Note that you define the z extent of the
   simulation box for 2d simulations when using the
   :doc:`create_box <create_box>` or :doc:`read_data <read_data>` commands.

**Output info:**

By default, this compute calculates a per-atom array with 2
columns. In regular dynamic tessellation mode the first column is the
Voronoi volume, the second is the neighbor count, as described above
(read above for the output data in case the *occupation* keyword is
specified).  These values can be accessed by any command that uses
per-atom values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options. If the *peratom* keyword is set to "no", the per-atom array
is still created, but it is not accessible.

If the *edge_histo* keyword is used, then this compute generates a
global vector of length *maxedge*\ +1, containing a histogram of the
number of edges per face.

If the *neighbors* value is set to yes, then this compute calculates a
local array with 3 columns. There is one row for each face of each
Voronoi cell.

.. note::

   Some LAMMPS commands such as the :doc:`compute reduce <compute_reduce>` command can accept either a per-atom or
   local quantity. If this compute produces both quantities, the command
   may access the per-atom quantity, even if you want to access the local
   quantity.  This effect can be eliminated by using the *peratom*
   keyword to turn off the production of the per-atom quantities.  For
   the default value *yes* both quantities are produced.  For the value
   *no*\ , only the local array is produced.

The Voronoi cell volume will be in distance :doc:`units <units>` cubed.
The Voronoi face area will be in distance :doc:`units <units>` squared.

Restrictions
""""""""""""

This compute is part of the VORONOI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

It also requires you have a copy of the Voro++ library built and
installed on your system.  See instructions on obtaining and
installing the Voro++ software in the src/VORONOI/README file.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`, :doc:`dump local <dump>`

**Default:** *neighbors* no, *peratom* yes
