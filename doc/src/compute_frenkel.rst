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
adding any atoms to the system).  Each atom in the compute group is then
assigned to the nearest lattice site.  Atoms that are not in the compute
group are ignored and never count as occupants of a lattice site; this
can be used, for example, to exclude gas atoms such as helium in
tungsten from the analysis.  A site with no atoms is a *vacancy*; a site
with two atoms holds an *interstitial*; a site with more than two atoms
is counted as an interstitial and additionally flagged as *irregular*.

Nearby defects are further grouped into clusters: a vacant site within
the distance *drvac*, or a multiply occupied site within the distance
*drint*, of another defective site belongs to the same cluster.  The
size of a cluster is the number of its interstitials minus the number of
its vacancies, so the sign of the size distinguishes interstitial-type
from vacancy-type clusters.  Note that the identification and counting
of the vacancies and interstitials themselves depends only on the number
of atoms at each site, not on this clustering or the two distances.

.. figure:: JPG/frenkel-diagram.png
   :figwidth: 50%
   :align: center

   Schematic depiction of a Frenkel pair: an atom is displaced
   from its lattice site leaving a vacancy and gets squeezed in
   between the atoms of neighboring occupied lattice sites

Several settings of this compute can be adjusted with the
:doc:`compute_modify <compute_modify>` command.  The following keywords
are recognized:

.. parsed-literal::

   *drvac* value = distance for including vacancies in a cluster (distance units)
   *drint* value = distance for including interstitials in a cluster (distance units)
   *region* value = ID of a region (or *none*) to restrict the lattice sites to
   *rescale* value = *yes* or *no* to co-scale the reference sites with the box
   *site_file* value = name of a file to read explicit "x y z" site coordinates from (or *none*)

The *drvac* and *drint* distances default to 1.01 and 1.42 lattice
spacings, respectively, so that vacancy clusters connect first- and
second-neighbor sites while the spatially more extended interstitial
configurations (e.g. dumbbells) connect over a somewhat larger distance.

Use the *region* keyword to exclude lattice sites where no atoms are
expected.  This is required when parts of the simulation box are empty,
for example the vacuum above a free surface in a non-periodic dimension;
otherwise all sites in the empty space are counted as vacancies.

Use *rescale yes* when the box changes size during the run (for example
under :doc:`fix npt <fix_nh>` or while heating), so that the reference
sites expand and contract with the box and thermal expansion is not
mistaken for defect formation.  A warning is printed when the box
changes during a run while *rescale* is not enabled.

With the *site_file* keyword the reference sites are not generated from
the lattice but read from the given text file, which must contain one
site per line as three coordinates "x y z" (anything following a "#"
character is ignored).  This allows using a reference structure that is
not a perfect lattice, for example the relaxed atom positions from the
beginning of the simulation.  A :doc:`lattice <lattice>` command is
still required, since the shape of the Wigner-Seitz cells and the
default distances are derived from it.

.. note::

   The lattice must match the crystal that the atoms actually occupy.  If
   the lattice spacing or orientation is wrong, essentially every atom
   will be flagged as a defect.  For a crystal at finite temperature it is
   usually best to use the thermally expanded lattice constant (or
   *rescale yes*), and to analyze the inherent structure (a quenched or
   minimized snapshot) when the thermal displacements are large.

In a restarted simulation this compute behaves like any other compute:
the :doc:`lattice <lattice>`, compute, and :doc:`compute_modify
<compute_modify>` commands must be repeated in the input script and the
reference sites are then regenerated at the beginning of the run.  With
*rescale yes* the lattice constant must be chosen to match the size of
the box stored in the restart file.

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
with a net content of *k* defects (clusters larger than 20 are added to
the last column).  Clusters containing as many vacancies as
interstitials, i.e. Frenkel pairs that are about to recombine, have a
net size of zero and are not included in the histogram.

The **per-atom vector** is the distance of each atom from its nearest
lattice site, which can be used to color atoms in a :doc:`dump image
<dump_image>` or to select displaced atoms.  For atoms outside the
compute group the value is 0.0.

The **local array** has 5 columns and one row per defect cluster owned by
the MPI process: the cluster ID, the cluster size (negative for vacancy
clusters, positive for interstitial clusters), and the *x*, *y*, *z*
coordinates of the cluster center.  To avoid double counting, a cluster
is stored only on the process whose subdomain contains its center.  The
array can be written with the :doc:`dump local <dump>` command.

The following excerpt from a displacement-cascade simulation in bcc iron
started by giving a 2 keV recoil to a primary knock-on atom (PKA) uses
compute frenkel to count and visualize the created Frenkel pairs.

.. code-block:: LAMMPS

   lattice      bcc 2.8553                     # reference lattice for the WS analysis
   compute      ke all ke/atom
   compute      fr all frenkel
   variable     vizstep index 100
   variable     hot atom c_ke>0.5
   variable     acol atom log(c_ke)
   group        hot dynamic all var hot every ${vizstep}

   fix          spheres all graphics/objects ${vizstep} &
                   sphere 1 32.0 10.0 0.0 2.0 &
                   sphere 2 50.0 20.0 0.0 2.0
   fix          label   all graphics/labels  ${vizstep} &
                   text "Frenkel pairs: $(c_fr[1]:% 4.0f)    Simulation time: $(time:% 5.1f) ps" &
                   400 50 0 size 30 &
                colorscale viz "log(kinetic energy / eV)" 700 400 0 vertical length 600 tics 10 &
                text "Vacancy" 280 185 0 size 30  &
                text "Interstitial" 300 115 0 size 30
   variable     acol atom log(c_ke)*v_hot

   dump         viz hot image 100 frenkel-*.png v_acol type size 800 800 &
                    zoom 2.0 view 70 30 center s 0.5 0.5 0.4 &
                    shiny 0.2 fsaa yes box no 0.0 axes yes 0.5 0.05 &
                    compute fr type 0 2 fix label type 1 0 fix spheres type 0 0

   dump_modify  viz pad 6 backcolor black backcolor2 white element Fe Fe Fe O &
                adiam * 2.5  color map2 0.342 0.062 0.429 color map3 0.736 0.216 0.330 &
                amap -0.5 1 cf 0.0 3 min map2 0.5 map3 max pink &
                acolor 1 steelblue acolor 2 darkgoldenrod

.. |frenkel1| image:: JPG/frenkel-sim-0.2.png
   :width: 33%

.. |frenkel2| image:: JPG/frenkel-sim-1.0.png
   :width: 33%

.. |frenkel3| image:: JPG/frenkel-sim-2.5.png
   :width: 33%


|frenkel1|  |frenkel2|  |frenkel3|

The images above are three snapshot images created by the LAMMPS input
from above.  Shown are atoms with elevated kinetic energy (smaller
spheres, colored by their kinetic energy on a logarithmic scale) and the
Frenkel pairs (larger spheres, blue: vacancies, yellow: interstitials).
The cascade of collisions spreads and briefly "melts" a small region
(0.2 ps) whose kinetic energy then quickly dissipates into the
surrounding crystal (1.0 ps) and the system relaxes and the lattice
reconstructs so that only a small number of Frenkel pairs survive.

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

Each cluster sphere is drawn with a diameter of 0.6 lattice spacings.
The *cflag2* setting is added to that diameter, which allows to enlarge
the markers; the *cflag1* setting is not used for spheres.


Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

All atoms must have IDs and an :doc:`atom map <atom_modify>` must be
defined (for example with ``atom_modify map array``).  A :doc:`lattice
<lattice>` must be defined to provide the reference sites; a general
triclinic lattice is not supported.

This compute cannot be used together with :doc:`comm_style tiled
<comm_style>` or :doc:`fix balance <fix_balance>`, and the *rescale*
option does not support triclinic simulation boxes.

Related commands
""""""""""""""""

:doc:`dump local <dump>`, :doc:`dump image <dump_image>`,
:doc:`compute cluster/atom <compute_cluster_atom>`,
:doc:`compute voronoi/atom <compute_voronoi_atom>`,
:doc:`lattice <lattice>`

Default
"""""""

The *drvac* and *drint* distances default to 1.01 and 1.42 lattice
spacings, respectively; *region* = none, *rescale* = no, *site_file* =
none.

----------

.. _compute-frenkel-Hammond:

**(Hammond)** Hammond, "Parallel point defect identification in molecular
dynamics simulations without post-processing: A compute and dump style for
LAMMPS", Comput. Phys. Commun. 247, 106862 (2020).
