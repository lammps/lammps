Use chunks to calculate system properties
=========================================

In LAMMS, "chunks" are collections of atoms, as defined by the
:doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns
each atom to a chunk ID (or to no chunk at all).  The number of chunks
and the assignment of chunk IDs to atoms can be static or change over
time.  Examples of "chunks" are molecules or spatial bins or atoms
with similar values (e.g. coordination number or potential energy).

The per-atom chunk IDs can be used as input to two other kinds of
commands, to calculate various properties of a system:

* :doc:`fix ave/chunk <fix_ave_chunk>`
* any of the :doc:`compute \*/chunk <compute>` commands

Here a brief overview for each of the 4 kinds of chunk-related commands
is provided.  Then some examples are given of how to compute different
properties with chunk commands.

Compute chunk/atom command:
---------------------------

This compute can assign atoms to chunks of various styles.  Only atoms
in the specified group and optional specified region are assigned to a
chunk.  Here are some possible chunk definitions:

+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| atoms in same molecule                                  | chunk ID = molecule ID                                                                                                          |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| atoms of same atom type                                 | chunk ID = atom type                                                                                                            |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| all atoms with same atom property (charge, radius, etc) | chunk ID = output of compute property/atom                                                                                      |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| atoms in same cluster                                   | chunk ID = output of :doc:`compute cluster/atom <compute_cluster_atom>` command                                                 |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| atoms in same spatial bin                               | chunk ID = bin ID                                                                                                               |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| atoms in same rigid body                                | chunk ID = molecule ID used to define rigid bodies                                                                              |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| atoms with similar potential energy                     | chunk ID = output of :doc:`compute pe/atom <compute_pe_atom>`                                                                   |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| atoms with same local defect structure                  | chunk ID = output of :doc:`compute centro/atom <compute_centro_atom>` or :doc:`compute coord/atom <compute_coord_atom>` command |
+---------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+

Note that chunk IDs are integer values, so for atom properties or
computes that produce a floating point value, they will be truncated
to an integer.  You could also use the compute in a variable that
scales the floating point value to spread it across multiple integers.

Spatial bins can be of various kinds, e.g. 1d bins = slabs, 2d bins =
pencils, 3d bins = boxes, spherical bins, cylindrical bins.

This compute also calculates the number of chunks *Nchunk*, which is
used by other commands to tally per-chunk data.  *Nchunk* can be a
static value or change over time (e.g. the number of clusters).  The
chunk ID for an individual atom can also be static (e.g. a molecule
ID), or dynamic (e.g. what spatial bin an atom is in as it moves).

Note that this compute allows the per-atom output of other
:doc:`computes <compute>`, :doc:`fixes <fix>`, and
:doc:`variables <variable>` to be used to define chunk IDs for each
atom.  This means you can write your own compute or fix to output a
per-atom quantity to use as chunk ID.  See the :doc:`Modify <Modify>`
doc pages for info on how to do this.  You can also define a :doc:`per-atom variable <variable>` in the input script that uses a formula to
generate a chunk ID for each atom.

Fix ave/chunk command:
----------------------

This fix takes the ID of a :doc:`compute chunk/atom <compute_chunk_atom>` command as input.  For each chunk,
it then sums one or more specified per-atom values over the atoms in
each chunk.  The per-atom values can be any atom property, such as
velocity, force, charge, potential energy, kinetic energy, stress,
etc.  Additional keywords are defined for per-chunk properties like
density and temperature.  More generally any per-atom value generated
by other :doc:`computes <compute>`, :doc:`fixes <fix>`, and :doc:`per-atom variables <variable>`, can be summed over atoms in each chunk.

Similar to other averaging fixes, this fix allows the summed per-chunk
values to be time-averaged in various ways, and output to a file.  The
fix produces a global array as output with one row of values per
chunk.

Compute \*/chunk commands:
--------------------------

The following computes operate on chunks of atoms to produce per-chunk
values.  Any compute whose style name ends in "/chunk" is in this
category:

* :doc:`compute com/chunk <compute_com_chunk>`
* :doc:`compute gyration/chunk <compute_gyration_chunk>`
* :doc:`compute inertia/chunk <compute_inertia_chunk>`
* :doc:`compute msd/chunk <compute_msd_chunk>`
* :doc:`compute property/chunk <compute_property_chunk>`
* :doc:`compute temp/chunk <compute_temp_chunk>`
* :doc:`compute torque/chunk <compute_vcm_chunk>`
* :doc:`compute vcm/chunk <compute_vcm_chunk>`

They each take the ID of a :doc:`compute chunk/atom <compute_chunk_atom>` command as input.  As their names
indicate, they calculate the center-of-mass, radius of gyration,
moments of inertia, mean-squared displacement, temperature, torque,
and velocity of center-of-mass for each chunk of atoms.  The :doc:`compute property/chunk <compute_property_chunk>` command can tally the
count of atoms in each chunk and extract other per-chunk properties.

The reason these various calculations are not part of the :doc:`fix ave/chunk command <fix_ave_chunk>`, is that each requires a more
complicated operation than simply summing and averaging over per-atom
values in each chunk.  For example, many of them require calculation
of a center of mass, which requires summing mass\*position over the
atoms and then dividing by summed mass.

All of these computes produce a global vector or global array as
output, with one or more values per chunk.  The output can be used in
various ways:

* As input to the :doc:`fix ave/time <fix_ave_time>` command, which can
  write the values to a file and optionally time average them.
* As input to the :doc:`fix ave/histo <fix_ave_histo>` command to
  histogram values across chunks.  E.g. a histogram of cluster sizes or
  molecule diffusion rates.
* As input to special functions of :doc:`equal-style variables <variable>`, like sum() and max() and ave().  E.g. to
  find the largest cluster or fastest diffusing molecule or average
  radius-of-gyration of a set of molecules (chunks).

Other chunk commands:
---------------------

* :doc:`compute chunk/spread/atom <compute_chunk_spread_atom>`
* :doc:`compute reduce/chunk <compute_reduce_chunk>`

The :doc:`compute chunk/spread/atom <compute_chunk_spread_atom>` command
spreads per-chunk values to each atom in the chunk, producing per-atom
values as its output.  This can be useful for outputting per-chunk
values to a per-atom :doc:`dump file <dump>`.  Or for using an atom's
associated chunk value in an :doc:`atom-style variable <variable>`.  Or
as input to the :doc:`fix ave/chunk <fix_ave_chunk>` command to
spatially average per-chunk values calculated by a per-chunk compute.

The :doc:`compute reduce/chunk <compute_reduce_chunk>` command reduces a
peratom value across the atoms in each chunk to produce a value per
chunk.  When used with the :doc:`compute chunk/spread/atom <compute_chunk_spread_atom>` command it can
create peratom values that induce a new set of chunks with a second
:doc:`compute chunk/atom <compute_chunk_atom>` command.

Example calculations with chunks
--------------------------------

Here are examples using chunk commands to calculate various
properties:

(1) Average velocity in each of 1000 2d spatial bins:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom bin/2d x 0.0 0.1 y lower 0.01 units reduced
   fix 1 all ave/chunk 100 10 1000 cc1 vx vy file tmp.out

(2) Temperature in each spatial bin, after subtracting a flow
velocity:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom bin/2d x 0.0 0.1 y lower 0.1 units reduced
   compute vbias all temp/profile 1 0 0 y 10
   fix 1 all ave/chunk 100 10 1000 cc1 temp bias vbias file tmp.out

(3) Center of mass of each molecule:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   compute myChunk all com/chunk cc1
   fix 1 all ave/time 100 1 100 c_myChunk[*] file tmp.out mode vector

(4) Total force on each molecule and ave/max across all molecules:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   fix 1 all ave/chunk 1000 1 1000 cc1 fx fy fz file tmp.out
   variable xave equal ave(f_1[2])
   variable xmax equal max(f_1[2])
   thermo 1000
   thermo_style custom step temp v_xave v_xmax

(5) Histogram of cluster sizes:

.. code-block:: LAMMPS

   compute cluster all cluster/atom 1.0
   compute cc1 all chunk/atom c_cluster compress yes
   compute size all property/chunk cc1 count
   fix 1 all ave/histo 100 1 100 0 20 20 c_size mode vector ave running beyond ignore file tmp.histo

(6) An example for using a per-chunk value to apply per-atom forces to
compress individual polymer chains (molecules) in a mixture, is
explained on the :doc:`compute chunk/spread/atom <compute_chunk_spread_atom>` command doc page.

(7) An example for using one set of per-chunk values for molecule
chunks, to create a second set of micelle-scale chunks (clustered
molecules, due to hydrophobicity), is explained on the
:doc:`compute reduce/chunk <compute_reduce_chunk>` command doc page.

(8) An example for using one set of per-chunk values (dipole moment
vectors) for molecule chunks, spreading the values to each atom in
each chunk, then defining a second set of chunks as spatial bins, and
using the :doc:`fix ave/chunk <fix_ave_chunk>` command to calculate an
average dipole moment vector for each bin.  This example is explained
on the :doc:`compute chunk/spread/atom <compute_chunk_spread_atom>`
command doc page.
