.. index:: create_box

create_box command
==================

Syntax
""""""

.. code-block:: LAMMPS

   create_box N region-ID keyword value ...
   create_box N NULL alo ahi blo bhi clo chi keyword value ...

* N = # of atom types to use in this simulation
* region-ID = ID of region to use as simulation domain or NULL for general triclinic box
* alo,ahi,blo,bhi,clo,chi = multipliers on a1,a2,a3 vectors defined by :doc"`lattice <lattice>` command (only when region-ID = NULL)
* zero or more keyword/value pairs may be appended
* keyword = *bond/types* or *angle/types* or *dihedral/types* or *improper/types* or *extra/bond/per/atom* or *extra/angle/per/atom* or *extra/dihedral/per/atom* or *extra/improper/per/atom* or *extra/special/per/atom*

  .. parsed-literal::

       *bond/types* value = # of bond types
       *angle/types* value = # of angle types
       *dihedral/types* value = # of dihedral types
       *improper/types* value = # of improper types
       *extra/bond/per/atom* value = # of bonds per atom
       *extra/angle/per/atom* value = # of angles per atom
       *extra/dihedral/per/atom* value = # of dihedrals per atom
       *extra/improper/per/atom* value = # of impropers per atom
       *extra/special/per/atom* value = # of special neighbors per atom

Examples
""""""""

.. code-block:: LAMMPS

   # orthogonal or restricted triclinic box using regionID = mybox
   create_box 2 mybox
   create_box 2 mybox bond/types 2 extra/bond/per/atom 1

.. code-block:: LAMMPS

   # 2d general triclinic box using primitive cell for 2d hex lattice
   lattice       custom 1.0 a1 1.0 0.0 0.0 a2 0.5 0.86602540378 0.0 &
                 a3 0.0 0.0 1.0 basis 0.0 0.0 0.0 triclinic/general
   create_box    1 NULL 0 5 0 5 -0.5 0.5

.. code-block:: LAMMPS

   # 3d general triclinic box using primitive cell for 3d fcc lattice
   lattice custom 1.0 a2 0.0 0.5 0.5 a1 0.5 0.0 0.5 a3 0.5 0.5 0.0 basis 0.0 0.0 0.0 triclinic/general
   create box 1 NULL -5 5 -10 10 0 20

Description
"""""""""""

This command creates a simulation box. It also partitions the box into
a regular 3d grid of smaller sub-boxes, one per processor (MPI task).
The geometry of the partitioning is based on the size and shape of the
simulation box, the number of processors being used and the settings
of the :doc:`processors <processors>` command.  The partitioning can
later be changed by the :doc:`balance <balance>` or :doc:`fix balance
<fix_balance>` commands.

Simulation boxes in LAMMPS can be either orthogonal or triclinic in
shape.  Orthogonal boxes are a brick in 3d (rectangle in 2d) with 6
faces that are each perpendicular to one of the standard xyz
coordinate axes.  Triclinic boxes are a parallelepiped in 3d
(parallelogram in 2d) with opposite pairs of faces parallel to each
other.  LAMMPS supports two forms of triclinic boxes, restricted and
general, which differ in how the box is oriented with respect to the
xyz coordinate axes.  See the :doc:`Howto triclinic <Howto_triclinic>`
for a detailed description of all 3 kinds of simulation boxes.

The argument *N* is the number of atom types that will be used in the
simulation.

Orthogonal and restricted triclinic boxes are created by specifying a
region ID previously defined by the :doc:`region <region>` command.
General triclinic boxes are discussed below.

If the region is not of style *prism*, then LAMMPS encloses the region
(block, sphere, etc.) with an axis-aligned orthogonal bounding box
which becomes the simulation domain.  For a 2d simulation, the zlo and
zhi values of the simulation box must straddle zero.

If the region is of style *prism*, LAMMPS creates a non-orthogonal
simulation domain shaped as a parallelepiped with triclinic symmetry.
As defined by the :doc:`region prism <region>` command, the
parallelepiped has an "origin" at (xlo,ylo,zlo) and three edge vectors
starting from the origin given by :math:`\vec a =
(x_\text{hi}-x_\text{lo},0,0)`; :math:`\vec b =
(xy,y_\text{hi}-y_\text{lo},0)`; and :math:`\vec c =
(xz,yz,z_\text{hi}-z_\text{lo})`.  In LAMMPS lingo, this is a
restricted triclinic box because the three edge vectors cannot be
defined in arbitrary (general) directions.  The parameters *xy*\ ,
*xz*\ , and *yz* can be 0.0 or positive or negative values and are
called "tilt factors" because they are the amount of displacement
applied to faces of an originally orthogonal box to transform it into
the parallelepiped.  For a 2d simulation, the zlo and zhi values of
the simulation box must straddle zero.

Typically a *prism* region used with the create_box command should
have tilt factors :math:`(xy,xz,yz)` that do not skew the box more
than half the distance of the parallel box length.  For example, if
:math:`x_\text{lo} = 2` and :math:`x_\text{hi} = 12`, then the
:math:`x` box length is 10 and the :math:`xy` tilt factor must be
between :math:`-5` and :math:`5`.  Similarly, both :math:`xz` and
:math:`yz` must be between :math:`-(x_\text{hi}-x_\text{lo})/2` and
:math:`+(y_\text{hi}-y_\text{lo})/2`.  Note that this is not a
limitation, since if the maximum tilt factor is 5 (as in this
example), then configurations with tilt :math:`= \dots, -15`,
:math:`-5`, :math:`5`, :math:`15`, :math:`25, \dots` are all
geometrically equivalent.

LAMMPS will issue a warning if the tilt factors of the created box do
not meet this criterion.  This is because simulations with large tilt
factors may run inefficiently, since they require more ghost atoms and
thus more communication.  With very large tilt factors, LAMMPS may
eventually produce incorrect trajectories and stop with errors due to
lost atoms or similar issues.

See the :doc:`Howto triclinic <Howto_triclinic>` page for geometric
descriptions of triclinic boxes and tilt factors, as well as how to
transform the restricted triclinic parameters to and from other
commonly used triclinic representations.

When a prism region is used, the simulation domain should normally be
periodic in the dimension that the tilt is applied to, which is given
by the second dimension of the tilt factor (e.g., :math:`y` for
:math:`xy` tilt).  This is so that pairs of atoms interacting across
that boundary will have one of them shifted by the tilt factor.
Periodicity is set by the :doc:`boundary <boundary>` command.  For
example, if the :math:`xy` tilt factor is non-zero, then the :math:`y`
dimension should be periodic.  Similarly, the :math:`z` dimension
should be periodic if :math:`xz` or :math:`yz` is non-zero.  LAMMPS
does not require this periodicity, but you may lose atoms if this is
not the case.

Note that if your simulation will tilt the box (e.g., via the
:doc:`fix deform <fix_deform>` command), the simulation box must be
created as triclinic, even if the tilt factors are initially 0.0.  You
can also change an orthogonal box to a triclinic box or vice versa by
using the :doc:`change box <change_box>` command with its *ortho* and
*triclinic* options.

.. note::

   If the system is non-periodic (in a dimension), then you should not
   make the lo/hi box dimensions (as defined in your :doc:`region
   <region>` command) radically smaller/larger than the extent of the
   atoms you eventually plan to create (e.g., via the
   :doc:`create_atoms <create_atoms>` command).  For example, if your
   atoms extend from 0 to 50, you should not specify the box bounds as
   :math:`-10000` and :math:`10000`. This is because as described
   above, LAMMPS uses the specified box size to lay out the 3d grid of
   processors.  A huge (mostly empty) box will be sub-optimal for
   performance when using "fixed" boundary conditions (see the
   :doc:`boundary <boundary>` command).  When using "shrink-wrap"
   boundary conditions (see the :doc:`boundary <boundary>` command), a
   huge (mostly empty) box may cause a parallel simulation to lose
   atoms the first time that LAMMPS shrink-wraps the box around the
   atoms.

----------

As noted above, general triclinic boxes in LAMMPS allow the box to
have arbitrary edge vectors **A**, **B**, **C**.  The only
restrictions are that the three vectors be distinct, non-zero, and not
co-planar.  They must also define a right-handed system such that
(**A** x **B**) points in the direction of **C**.  Note that a
left-handed system can be converted to a right-handed system by simply
swapping the order of any pair of the **A**, **B**, **C** vectors.

To create a general triclinic boxes, the region is specified as NULL
and the next 6 parameters (alo,ahi,blo,bhi,clo,chi) define the three
edge vectors **A**, **B**, **C** using additional information
previously defined by the :doc:`lattice <lattice>` command.

The lattice must be of style *custom* and use its *triclinic/general*
option.  This insures the lattice satisfies the restrictions listed
above.  The *a1, *a2*, *a3* settings of the :doc:`lattice <lattice>`
command define the edge vectors of a unit cell of the general
triclinic lattice.  This command uses them to define the three edge
vectors and origin of the general triclinic box as:

* **A** = (ahi-alo) * *a1*
* **B** = (bhi-blo) * *a2*
* **C** = (chi-clo) * *a3*
* origin = (alo*a1 + blo*a2 + clo*a3)

For 2d general triclinic boxes, clo = -0.5 and chi = 0.5 is required.

.. note::

   LAMMPS allows specification of general triclinic simulation boxes
   as a convenience for users who may be converting data from
   solid-state crystallographic representations or from DFT codes for
   input to LAMMPS.  However, as explained on the
   :doc:`Howto_triclinic <Howto_triclinic>` doc page, internally,
   LAMMPS only uses restricted triclinic simulation boxes.  This means
   the box defined by this command and per-atom information
   (e.g. coordinates, velocities) defined by the :doc:`create_atoms
   <create_atoms>` command are converted (rotated) from general to
   restricted triclinic form when the two commands are invoked.  The
   <Howto_triclinic>` doc page also discusses other LAMMPS commands
   which can input/output general triclinic representations of the
   simulation box and per-atom data.

----------

The optional keywords can be used to create a system that allows for
bond (angle, dihedral, improper) interactions, or for molecules with
special 1--2, 1--3, or 1--4 neighbors to be added later.  These
optional keywords serve the same purpose as the analogous keywords
that can be used in a data file which are recognized by the
:doc:`read_data <read_data>` command when it sets up a system.

Note that if these keywords are not used, then the create_box command
creates an atomic (non-molecular) simulation that does not allow bonds
between pairs of atoms to be defined, or a :doc:`bond potential
<bond_style>` to be specified, or for molecules with special neighbors
to be added to the system by commands such as :doc:`create_atoms mol
<create_atoms>`, :doc:`fix deposit <fix_deposit>` or :doc:`fix pour
<fix_pour>`.

As an example, see the examples/deposit/in.deposit.molecule script,
which deposits molecules onto a substrate.  Initially there are no
molecules in the system, but they are added later by the :doc:`fix
deposit <fix_deposit>` command.  The create_box command in the script
uses the bond/types and extra/bond/per/atom keywords to allow this.
If the added molecule contained more than one special bond (allowed by
default), an extra/special/per/atom keyword would also need to be
specified.

----------

Restrictions
""""""""""""

An :doc:`atom_style <atom_style>` and :doc:`region <region>` must have
been previously defined to use this command.

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`create_atoms <create_atoms>`,
:doc:`region <region>`

Default
"""""""

none
