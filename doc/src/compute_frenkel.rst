.. index:: compute frenkel

compute frenkel command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID frenkel

* ID, group-ID are documented in :doc:`compute <compute>` command
* frenkel = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all frenkel

   compute def all frenkel
   compute_modify def drvac 1.2 drint 1.7 rescale yes

Description
"""""""""""

.. versionadded:: TBD

Define a computation that identifies point defects in a crystal by
Wigner-Seitz analysis and counts the number of vacancies and
interstitials (i.e.\ the number of Frenkel pairs) on the fly, without
post-processing.  This is useful for following radiation-damage cascades,
thermal defect formation, or any process that creates or annihilates
point defects.

The reference lattice is taken from the most recently defined
:doc:`lattice <lattice>` command: the compute generates one lattice site
at every basis point of every unit cell that lies inside the simulation
box (exactly as :doc:`create_atoms <create_atoms>` would, but without
adding any atoms to the system).  Each atom is then assigned to the
nearest lattice site.  A site with no atoms is a *vacancy*; a site with
two atoms holds an *interstitial*; a site with more than two atoms is
counted as an interstitial and additionally flagged as *irregular*.
Neighboring vacant or over-occupied sites are grouped into clusters.

The search radii used to associate atoms with sites can be adjusted with
the *drvac* and *drint* keywords of the :doc:`compute_modify
<compute_modify>` command; by default they are derived from the lattice
spacing.  The following :doc:`compute_modify <compute_modify>` keywords
are recognized:

.. parsed-literal::

   *drvac* value = cutoff distance for detecting vacancies (distance units)
   *drint* value = cutoff distance for detecting interstitials (distance units)
   *region* value = ID of a region (or *none*) to restrict the lattice sites to
   *frenkelgroup* value = ID of a group of atoms to restrict the defect search to
   *rescale* value = *yes* or *no* to co-scale the reference sites with the box
   *site_file* value = name of a file with explicit "x y z" site coordinates (or *none*)

Use *rescale yes* when the box changes size during the run (for example
under :doc:`fix npt <fix_nh>` or while heating), so that the reference
sites expand and contract with the box and thermal expansion is not
mistaken for defect formation.

.. note::

   The lattice must match the crystal that the atoms actually occupy.  If
   the lattice spacing or orientation is wrong, essentially every atom
   will be flagged as a defect.  For a crystal at finite temperature it is
   usually best to use the thermally expanded lattice constant (or
   *rescale yes*), and to analyze the inherent structure (a quenched or
   minimized snapshot) when the thermal displacements are large.

This compute is described in :ref:`(Hammond) <compute-frenkel-Hammond>`.

Output info
"""""""""""

This compute calculates a global vector of length 3, a global array, a
per-atom vector, and a local array.

The **global vector** holds, in order, the number of vacancies, the
number of interstitials, and the number of irregular sites (sites with
more than two atoms), each summed over all MPI processes.  Thus
``c_ID[1]`` is the number of Frenkel pairs.

The **global array** has 2 rows and 20 columns and is a histogram of the
defect cluster sizes: row 1 counts vacancy clusters and row 2 counts
interstitial clusters, with column *k* holding the number of clusters
containing *k* defects (clusters larger than 20 are added to the last
column).

The **per-atom vector** is the distance of each atom from its nearest
lattice site, which can be used to color atoms in a :doc:`dump image
<dump_image>` or to select displaced atoms.

The **local array** has 5 columns and one row per defect cluster owned by
the MPI process: the cluster ID, the cluster size (negative for vacancy
clusters, positive for interstitial clusters), and the *x*, *y*, *z*
coordinates of the cluster center.  To avoid double counting, a cluster
is stored only on the process whose subdomain contains its center.  The
array can be written with the :doc:`dump local <dump>` command.

The following excerpt from a displacement-cascade simulation in bcc iron
(a primary knock-on atom has been given a high velocity in the steps
before) follows the defect count on the fly and records the surviving
defects without any post-processing:

.. code-block:: LAMMPS

   lattice         bcc 2.8553                     # reference lattice for the WS analysis
   compute         fr all frenkel
   # live vacancy / interstitial / irregular-site count in the thermo output
   thermo_style    custom step time temp c_fr[1] c_fr[2] c_fr[3]
   # log the number of Frenkel pairs c_fr[1] versus time
   variable        t equal time
   fix             rec all ave/time 20 1 20 v_t c_fr[1] file frenkel.dat
   # write the surviving defects (size, x, y, z) every 1000 steps
   dump            dl all local 1000 defects.dump index c_fr[2] c_fr[3] c_fr[4] c_fr[5]
   run             15000

Here ``c_fr[1]`` in the :doc:`thermo_style <thermo_style>` and :doc:`fix
ave/time <fix_ave_time>` commands refers to the *global* vector (the
defect counts), while ``c_fr[2]`` ... ``c_fr[5]`` in the :doc:`dump local
<dump>` command refer to the columns of the *local* array (defect size and
position).

Dump image info
"""""""""""""""

Compute *frenkel* can be used with the *compute* keyword of :doc:`dump
image <dump_image>`.  It adds one sphere at the center of every defect
cluster to the rendered image, so the spatial distribution of the damage
is shown directly without an external visualization tool.

Each sphere carries a color index of 1 for a vacancy cluster and 2 for an
interstitial cluster.  With color style *type* or *element* these indices
are mapped to the corresponding atom-type (or element) colors; with color
style *const* all spheres use one color, which defaults to white and can
be changed with :doc:`dump_modify ccolor <dump_image>`.  The opacity
defaults to fully opaque and can be changed with *dump_modify ctrans*.

To draw vacancies and interstitials in two distinct colors that are
independent of the real atoms, define one or two extra atom types that no
atoms actually use (give them a mass and a :doc:`pair_coeff <pair_coeff>`
so the input is valid; since no atoms have those types they have no effect
on the simulation) and color them with :doc:`dump_modify acolor
<dump_image>`.  For example, with the metal atoms on type 3 and types 1
and 2 reserved for the defect colors:

.. code-block:: LAMMPS

   compute fr all frenkel
   dump    d all image 1000 defect.*.jpg type type adiam 0.5 compute fr type 0 0
   dump_modify d acolor 1 blue acolor 2 red acolor 3 gray atrans 3 0.1

The *cflag1* setting is added to the sphere diameter (the compute does not
know a meaningful radius, so it is set to 0; use a positive *cflag2*).
The *cflag2* setting is currently unused for spheres.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

All atoms must have IDs and an :doc:`atom map <atom_modify>` must be
defined (for example with ``atom_modify map array``).  A :doc:`lattice
<lattice>` must be defined to provide the reference sites; a general
triclinic lattice is not supported.

Related commands
""""""""""""""""

:doc:`dump local <dump>`, :doc:`dump image <dump_image>`,
:doc:`compute cluster/atom <compute_cluster_atom>`,
:doc:`compute voronoi/atom <compute_voronoi_atom>`,
:doc:`lattice <lattice>`

Default
"""""""

The *drvac* and *drint* cutoffs default to values derived from the lattice
spacing; *region* = none, *frenkelgroup* = the compute group, *rescale* =
no, *site_file* = none.

----------

.. _compute-frenkel-Hammond:

**(Hammond)** Hammond, "Parallel point defect identification in molecular
dynamics simulations without post-processing: A compute and dump style for
LAMMPS", Comput. Phys. Commun. 247, 106862 (2020).
