.. index:: dump
.. index:: dump atom
.. index:: dump cfg
.. index:: dump custom
.. index:: dump dcd
.. index:: dump grid
.. index:: dump grid/vtk
.. index:: dump local
.. index:: dump xtc
.. index:: dump yaml
.. index:: dump xyz
.. index:: dump atom/gz
.. index:: dump cfg/gz
.. index:: dump custom/gz
.. index:: dump local/gz
.. index:: dump xyz/gz
.. index:: dump atom/zstd
.. index:: dump cfg/zstd
.. index:: dump custom/zstd
.. index:: dump xyz/zstd
.. index:: dump local/zstd

dump command
============

:doc:`dump vtk <dump_vtk>` command
==================================

:doc:`dump h5md <dump_h5md>` command
====================================

:doc:`dump molfile <dump_molfile>` command
==========================================

:doc:`dump netcdf <dump_netcdf>` command
========================================

:doc:`dump image <dump_image>` command
======================================

:doc:`dump movie <dump_image>` command
======================================

:doc:`dump atom/adios <dump_adios>` command
===========================================

:doc:`dump custom/adios <dump_adios>` command
=============================================

:doc:`dump cfg/uef <dump_cfg_uef>` command
==========================================

Syntax
""""""

.. code-block:: LAMMPS

   dump ID group-ID style N file attribute1 attribute2 ...

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* style = *atom* or *atom/adios* or *atom/gz* or *atom/zstd* or *cfg* or *cfg/gz* or *cfg/zstd* or *cfg/uef* or *custom* or *custom/gz* or *custom/zstd* or *custom/adios* or *dcd* or *grid* or *grid/vtk* or *h5md* or *image* or *local* or *local/gz* or *local/zstd* or *molfile* or *movie* or *netcdf* or *netcdf/mpiio* or *vtk* or *xtc* or *xyz* or *xyz/gz* or *xyz/zstd* or *yaml*
* N = dump on timesteps which are multiples of N
* file = name of file to write dump info to
* attribute1,attribute2,... = list of attributes for a particular style

  .. parsed-literal::

       *atom* attributes = none
       *atom/adios* attributes = none,  discussed on :doc:`dump atom/adios <dump_adios>` page
       *atom/gz* attributes = none
       *atom/zstd* attributes = none
       *cfg* attributes = same as *custom* attributes, see below
       *cfg/gz* attributes = same as *custom* attributes, see below
       *cfg/zstd* attributes = same as *custom* attributes, see below
       *cfg/uef* attributes = same as *custom* attributes, discussed on :doc:`dump cfg/uef <dump_cfg_uef>` page
       *custom*, *custom/gz*, *custom/zstd* attributes = see below
       *custom/adios* attributes = same as *custom* attributes, discussed on :doc:`dump custom/adios <dump_adios>` page
       *dcd* attributes = none
       *h5md* attributes = discussed on :doc:`dump h5md <dump_h5md>` page
       *grid* attributes = see below
       *grid/vtk* attributes = see below
       *image* attributes = discussed on :doc:`dump image <dump_image>` page
       *local*, *local/gz*, *local/zstd* attributes = see below
       *molfile* attributes = discussed on :doc:`dump molfile <dump_molfile>` page
       *movie* attributes = discussed on :doc:`dump image <dump_image>` page
       *netcdf* attributes = discussed on :doc:`dump netcdf <dump_netcdf>` page
       *netcdf/mpiio* attributes = discussed on :doc:`dump netcdf <dump_netcdf>` page
       *vtk* attributes = same as *custom* attributes, see below, also :doc:`dump vtk <dump_vtk>` page
       *xtc* attributes = none
       *xyz* attributes = none
       *xyz/gz* attributes = none
       *xyz/zstd* attributes = none
       *yaml* attributes = same as *custom* attributes, see below

* *custom* or *custom/gz* or *custom/zstd* or *cfg* or *cfg/gz* or *cfg/zstd* or *cfg/uef* or *netcdf* or *netcdf/mpiio* or *yaml* attributes:

  .. parsed-literal::

         possible attributes = id, mol, proc, procp1, type, element, mass,
                               x, y, z, xs, ys, zs, xu, yu, zu,
                               xsu, ysu, zsu, ix, iy, iz,
                               vx, vy, vz, fx, fy, fz,
                               q, mux, muy, muz, mu,
                               radius, diameter, omegax, omegay, omegaz,
                               angmomx, angmomy, angmomz, tqx, tqy, tqz,
                               c_ID, c_ID[I], f_ID, f_ID[I], v_name,
                               i_name, d_name, i2_name[I], d2_name[I]

  .. parsed-literal::

           id = atom ID
           mol = molecule ID
           proc = ID of processor that owns atom
           procp1 = ID+1 of processor that owns atom
           type = atom type
           element = name of atom element, as defined by :doc:`dump_modify <dump_modify>` command
           mass = atom mass
           x,y,z = unscaled atom coordinates
           xs,ys,zs = scaled atom coordinates
           xu,yu,zu = unwrapped atom coordinates
           xsu,ysu,zsu = scaled unwrapped atom coordinates
           ix,iy,iz = box image that the atom is in
           vx,vy,vz = atom velocities
           fx,fy,fz = forces on atoms
           q = atom charge
           mux,muy,muz = orientation of dipole moment of atom
           mu = magnitude of dipole moment of atom
           radius,diameter = radius, diameter of spherical particle
           omegax,omegay,omegaz = angular velocity of spherical particle
           angmomx,angmomy,angmomz = angular momentum of aspherical particle
           tqx,tqy,tqz = torque on finite-size particles
           c_ID = per-atom vector calculated by a compute with ID
           c_ID[I] = Ith column of per-atom array calculated by a compute with ID, I can include wildcard (see below)
           f_ID = per-atom vector calculated by a fix with ID
           f_ID[I] = Ith column of per-atom array calculated by a fix with ID, I can include wildcard (see below)
           v_name = per-atom vector calculated by an atom-style variable with name
           i_name = custom integer vector with name
           d_name = custom floating point vector with name
           i2_name[I] = Ith column of custom integer array with name, I can include wildcard (see below)
           d2_name[I] = Ith column of custom floating point vector with name, I can include wildcard (see below)

* *local* or *local/gz* or *local/zstd* attributes:

  .. parsed-literal::

         possible attributes = index, c_ID, c_ID[I], f_ID, f_ID[I]
           index = enumeration of local values
           c_ID = local vector calculated by a compute with ID
           c_ID[I] = Ith column of local array calculated by a compute with ID, I can include wildcard (see below)
           f_ID = local vector calculated by a fix with ID
           f_ID[I] = Ith column of local array calculated by a fix with ID, I can include wildcard (see below)

* *grid* or *grid/vtk* attributes:

  .. parsed-literal::

         possible attributes = c_ID:gname:dname, c_ID:gname:dname[I], f_ID:gname:dname, f_ID:gname:dname[I]
           gname = name of grid defined by compute or fix
           dname = name of data field defined by compute or fix
           c_ID = per-grid vector calculated by a compute with ID
           c_ID[I] = Ith column of per-grid array calculated by a compute with ID, I can include wildcard (see below)
           f_ID = per-grid vector calculated by a fix with ID
           f_ID[I] = Ith column of per-grid array calculated by a fix with ID, I can include wildcard (see below)

Examples
""""""""

.. code-block:: LAMMPS

   dump myDump all atom 100 dump.lammpstrj
   dump myDump all atom/gz 100 dump.atom.gz
   dump myDump all atom/zstd 100 dump.atom.zst
   dump 2 subgroup atom 50 dump.run.bin
   dump 4a all custom 100 dump.myforce.* id type x y vx fx
   dump 4a all custom 100 dump.myvel.lammpsbin id type x y z vx vy vz
   dump 4b flow custom 100 dump.%.myforce id type c_myF[3] v_ke
   dump 4b flow custom 100 dump.%.myforce id type c_myF[*] v_ke
   dump 2 inner cfg 10 dump.snap.*.cfg mass type xs ys zs vx vy vz
   dump snap all cfg 100 dump.config.*.cfg mass type xs ys zs id type c_Stress[2]
   dump 1 all xtc 1000 file.xtc

Description
"""""""""""

Dump a snapshot of quantities to one or more files once every
:math:`N` timesteps in one of several styles.  The timesteps on which
dump output is written can also be controlled by a variable.  See the
:doc:`dump_modify every <dump_modify>` command.

Almost all the styles output per-atom data, i.e. one or more values
per atom.  The exceptions are as follows.  The *local* styles output
one or more values per bond (angle, dihedral, improper) or per pair of
interacting atoms (force or neighbor interactions).  The *grid* styles
output one or more values per grid cell, which are produced by other
commands which overlay the simulation domain with a regular grid.  See
the :doc:`Howto grid <Howto_grid>` doc page for details.  The *image*
style renders a JPG, PNG, or PPM image file of the system for each
snapshot, while the *movie* style combines and compresses the series
of images into a movie file; both styles are discussed in detail on
the :doc:`dump image <dump_image>` page.

Only information for atoms in the specified group is dumped.  The
:doc:`dump_modify thresh and region and refresh <dump_modify>`
commands can also alter what atoms are included.  Not all styles
support these options; see details on the :doc:`dump_modify
<dump_modify>` doc page.

As described below, the filename determines the kind of output: text
or binary or gzipped, one big file or one per timestep, one file for
all the processors or multiple smaller files.

.. note::

   Because periodic boundary conditions are enforced only on timesteps
   when neighbor lists are rebuilt, the coordinates of an atom written
   to a dump file may be slightly outside the simulation box.
   Re-neighbor timesteps will not typically coincide with the
   timesteps dump snapshots are written.  See the :doc:`dump_modify
   pbc <dump_modify>` command if you wish to force coordinates to be
   strictly inside the simulation box.

.. note::

   Unless the :doc:`dump_modify sort <dump_modify>` option is invoked,
   the lines of atom or grid information written to dump files
   (typically one line per atom or grid cell) will be in an
   indeterminate order for each snapshot.  This is even true when
   running on a single processor, if the :doc:`atom_modify sort
   <atom_modify>` option is on, which it is by default.  In this case
   atoms are re-ordered periodically during a simulation, due to
   spatial sorting.  It is also true when running in parallel, because
   data for a single snapshot is collected from multiple processors,
   each of which owns a subset of the atoms.

.. warning::

   Without either including atom IDs or using the :doc:`dump_modify sort
   <dump_modify>` option, it is impossible for visualization programs
   (e.g. OVITO or VMD) or analysis tools to assign data in different
   frames consistently to the same atom.  This can lead to incorrect
   visualizations or results.  LAMMPS will print a warning in such cases.

For the *atom*, *custom*, *cfg*, *grid*, and *local* styles, sorting
is off by default.  For the *dcd*, *grid/vtk*, *xtc*, *xyz*, and
*molfile* styles, sorting by atom ID or grid ID is on by default. See
the :doc:`dump_modify <dump_modify>` page for details.

The *style* keyword determines what kind of data is written to the
dump file(s) and in what format.

Note that *atom*, *custom*, *dcd*, *xtc*, *xyz*, and *yaml* style dump
files can be read directly by `VMD <https://www.ks.uiuc.edu/Research/vmd>`_,
a popular tool for visualizing and analyzing trajectories from atomic
and molecular systems.  For reading *netcdf* style dump files, the
netcdf plugin needs to be recompiled from source using a NetCDF version
compatible with the one used by LAMMPS.  The bundled plugin binary
uses a very old version of NetCDF that is not compatible with LAMMPS.

Likewise the `OVITO visualization package <https://www.ovito.org>`_,
popular for materials modeling, can read the *atom*, *custom*,
*local*, *xtc*, *cfg*, *netcdf*, and *xyz* style atom dump files
directly.  With version 3.8 and above, OVITO can also read and
visualize *grid* style dump files with grid cell data, including
iso-surface images of the grid cell values.

Note that settings made via the :doc:`dump_modify <dump_modify>`
command can also alter the format of individual values and content of
the dump file itself.  This includes the precision of values output to
text-based dump files which is controlled by the :doc:`dump_modify
format <dump_modify>` command and its options.

----------

Format of native LAMMPS format dump files:

The *atom*, *custom*, *grid*, and *local* styles create files in a
simple LAMMPS-specific text format that is mostly self-explanatory
when viewing a dump file.  Many post-processing tools either included
with LAMMPS or third-party tools can read this format, as does the
:doc:`rerun <rerun>` command.  See tools described on the :doc:`Tools
<Tools>` doc page for examples, including `Pizza.py
<https://lammps.github.io/pizza>`_.

For all these styles, the dimensions of the simulation box are
included in each snapshot.  The simulation box in LAMMPS can be
defined in one of 3 ways: orthogonal, restricted triclinic, and
general triclinic.  See the :doc:`Howto triclinic <Howto_triclinic>`
doc page for a detailed description of all 3 options.

For an orthogonal simulation box the box information is formatted as:

.. parsed-literal::

   ITEM: BOX BOUNDS xx yy zz
   xlo xhi
   ylo yhi
   zlo zhi

where xlo,xhi are the maximum extents of the simulation box in the
:math:`x`-dimension, and similarly for :math:`y` and :math:`z`.  The
"xx yy zz" terms are six characters that encode the style of boundary for each
of the six simulation box boundaries (xlo,xhi; ylo,yhi; and zlo,zhi).  Each of
the six characters is one of *p* (periodic), *f* (fixed), *s* (shrink wrap),
or *m* (shrink wrapped with a minimum value).  See the
:doc:`boundary <boundary>` command for details.

For a restricted triclinic simulation box, an orthogonal bounding box
which encloses the restricted triclinic simulation box is output,
along with the three tilt factors (*xy*, *xz*, *yz*) of the triclinic
box, formatted as follows:

.. parsed-literal::

   ITEM: BOX BOUNDS xy xz yz xx yy zz
   xlo_bound xhi_bound xy
   ylo_bound yhi_bound xz
   zlo_bound zhi_bound yz

The presence of the text "xy xz yz" in the ITEM line indicates that
the three tilt factors will be included on each of the three following lines.
This bounding box is convenient for many visualization programs.  The
meaning of the six character flags for "xx yy zz" is the same as above.

Note that the first two numbers on each line are now xlo_bound instead
of xlo, etc. because they represent a bounding box.  See the :doc:`Howto
triclinic <Howto_triclinic>` page for a geometric description of
triclinic boxes, as defined by LAMMPS, simple formulas for how the six
bounding box extents (xlo_bound, xhi_bound, etc.) are calculated from the
triclinic parameters, and how to transform those parameters to and
from other commonly used triclinic representations.

For a general triclinic simulation box, see the "General triclinic"
section below for a description of the ITEM: BOX BOUNDS format as well
as how per-atom coordinates and per-atom vector quantities are output.

The *atom* and *custom* styles output a "ITEM: NUMBER OF ATOMS" line
with the count of atoms in the snapshot.  Likewise they output an
"ITEM: ATOMS" line which includes column descriptors for the per-atom
lines that follow.  For example, the descriptors would be "id type xs
ys zs" for the default *atom* style, and would be the atom attributes
you specify in the dump command for the *custom* style.  Each
subsequent line will list the data for a single atom.

For style *atom*, atom coordinates are written to the file, along with
the atom ID and atom type.  By default, atom coords are written in a
scaled format (from 0 to 1).  That is, an :math:`x` value of 0.25 means the
atom is at a location 1/4 of the distance from *xlo* to *xhi* of the box
boundaries.  The format can be changed to unscaled coords via the
:doc:`dump_modify <dump_modify>` settings.  Image flags can also be
added for each atom via dump_modify.

Style *custom* allows you to specify a list of atom attributes to be
written to the dump file for each atom.  Possible attributes are
listed above and will appear in the order specified.  You cannot
specify a quantity that is not defined for a particular
simulation---such as *q* for atom style *bond*, since that atom style
does not assign charges.  Dumps occur at the very end of a timestep,
so atom attributes will include effects due to fixes that are applied
during the timestep.  An explanation of the possible dump custom
attributes is given below.

.. versionadded:: 22Dec2022

For style *grid* the dimension of the simulation domain and size of
the Nx by Ny by Nz grid that overlays the simulation domain are also
output with each snapshot:

.. parsed-literal::

   ITEM: DIMENSION
   dim
   ITEM: GRID SIZE
   nx ny nz

The value dim will be 2 or 3 for 2d or 3d simulations.  It is included
so that post-processing tools like `OVITO <https://www.ovito.org>`_,
which can visualize grid-based quantities know how to draw each grid
cell.  The grid size will match the input script parameters for
grid(s) created by the computes or fixes which are referenced by the
the dump command.  For 2d simulations (and grids), nz will always
be 1.

There will also be an "ITEM: GRID DATA" line which includes column
descriptors for the per grid cell data.  Each subsequent line (Nx *
Ny * Nz lines) will list the data for a single grid cell.  If grid
cell IDs are included in the output via the :doc:`compute
property/grid <compute_property_grid>` command, then the IDs will
range from 1 to N = Nx*Ny*Nz.  The ordering of IDs is with the x index
varying fastest, then the y index, and the z index varying slowest.

For style *local*, local output generated by :doc:`computes <compute>`
and :doc:`fixes <fix>` is used to generate lines of output that is
written to the dump file.  This local data is typically calculated by
each processor based on the atoms it owns, but there may be zero or
more entities per atom (e.g., a list of bond distances).  An explanation
of the possible dump local attributes is given below.  Note that by
using input from the :doc:`compute property/local
<compute_property_local>` command with dump local, it is possible to
generate information on bonds, angles, etc. that can be cut and pasted
directly into a data file read by the :doc:`read_data <read_data>`
command.

----------

Dump files in other popular formats:

.. note::

   This section only discusses file formats relevant to this doc page.
   The top of this page has links to other dump commands (with their
   own pages) which write files in additional popular formats.

Style *cfg* has the same command syntax as style *custom* and writes
extended CFG format files, as used by the `AtomEye
<http://li.mit.edu/Archive/Graphics/A/>`_ visualization package.
Since the extended CFG format uses a single snapshot of the system per
file, a wildcard "\*" must be included in the filename, as discussed
below.  The list of atom attributes for style *cfg* must begin with
either "mass type xs ys zs" or "mass type xsu ysu zsu" since these
quantities are needed to write the CFG files in the appropriate format
(though the "mass" and "type" fields do not appear explicitly in the
file).  Any remaining attributes will be stored as "auxiliary
properties" in the CFG files.  Note that you will typically want to
use the :doc:`dump_modify element <dump_modify>` command with
CFG-formatted files, to associate element names with atom types, so
that AtomEye can render atoms appropriately. When unwrapped
coordinates *xsu*, *ysu*, and *zsu* are requested, the nominal AtomEye
periodic cell dimensions are expanded by a large factor UNWRAPEXPAND =
10.0, which ensures atoms that are displayed correctly for up to
UNWRAPEXPAND/2 periodic boundary crossings in any direction.  Beyond
this, AtomEye will rewrap the unwrapped coordinates.  The expansion
causes the atoms to be drawn farther away from the viewer, but it is
easy to zoom the atoms closer, and the interatomic distances are
unaffected.

The *dcd* style writes DCD files, a standard atomic trajectory format
used by the CHARMM, NAMD, and XPlor molecular dynamics packages.  DCD
files are binary and thus may not be portable to different machines.
The number of atoms per snapshot cannot change with the *dcd* style.
The *unwrap* option of the :doc:`dump_modify <dump_modify>` command
allows DCD coordinates to be written "unwrapped" by the image flags
for each atom.  Unwrapped means that if the atom has passed through
a periodic boundary one or more times, the value is printed for what
the coordinate would be if it had not been wrapped back into the
periodic box.  Note that these coordinates may thus be far outside
the box size stored with the snapshot.

The *xtc* style writes XTC files, a compressed trajectory format used
by the GROMACS molecular dynamics package, and described `here
<https://manual.gromacs.org/current/reference-manual/file-formats.html#xtc>`_.
The precision used in XTC files can be adjusted via the
:doc:`dump_modify <dump_modify>` command.  The default value of 1000
means that coordinates are stored to 1/1000 nanometer accuracy.  XTC
files are portable binary files written in the NFS XDR data format, so
that any machine which supports XDR should be able to read them.  The
number of atoms per snapshot cannot change with the *xtc* style.  The
*unwrap* option of the :doc:`dump_modify <dump_modify>` command allows
XTC coordinates to be written "unwrapped" by the image flags for each
atom.  Unwrapped means that if the atom has passed through a periodic
boundary one or more times, the value is printed for what the
coordinate would be if it had not been wrapped back into the periodic
box.  Note that these coordinates may thus be far outside the box size
stored with the snapshot.

The *xyz* style writes XYZ files, which is a simple text-based
coordinate format that many codes can read. Specifically it has a line
with the number of atoms, then a comment line that is usually ignored
followed by one line per atom with the atom type and the :math:`x`-,
:math:`y`-, and :math:`z`-coordinate of that atom.  You can use the
:doc:`dump_modify element <dump_modify>` option to change the output
from using the (numerical) atom type to an element name (or some other
label). This will help many visualization programs to guess bonds and
colors.

.. versionadded:: 22Dec2022

The *grid/vtk* style writes VTK files for grid data on a regular
rectilinear grid.  Its content is conceptually similar to that of the
text file produced by the *grid* style, except that it in an XML-based
format which visualization programs which support the VTK format can
read, e.g. the `ParaView tool <https://www.paraview.org>`_.  For this
style, there can only be 1 or 3 per grid cell attributes specified.
If it is a single value, it is a scalar quantity.  If 3 values are
specified it is encoded in the VTK file as a vector quantity (for each
grid cell).  The filename for this style must include a "\*" wildcard
character to produce one file per snapshot; see details below.

.. versionadded:: 4May2022

Dump style *yaml* has the same command syntax as style *custom* and
writes YAML format files that can be easily parsed by a variety of
data processing tools and programming languages.  Each timestep will
be written as a YAML "document" (i.e., starts with "---" and ends with
"...").  The style supports writing one file per timestep through the
"\*" wildcard but not multi-processor outputs with the "%" token in
the filename.  In addition to per-atom data, :doc:`thermo <thermo>`
data can be included in the *yaml* style dump file using the
:doc:`dump_modify thermo yes <dump_modify>`. The data included in the
dump file uses the "thermo" tag and is otherwise identical to data
specified by the :doc:`thermo_style <thermo_style>` command.

Below is an example for a YAML format dump created by the following commands.

.. code-block:: LAMMPS

   dump out all yaml 100 dump.yaml id type x y z vx vy vz ix iy iz
   dump_modify out time yes units yes thermo yes format 1 %5d format "% 10.6e"

The tags "time", "units", and "thermo" are optional and enabled by the
dump_modify command. The list under the "box" tag has three lines for
orthogonal boxes and four lines for triclinic boxes, where the first
three are the box boundaries and the fourth the three tilt factors
(:math:`xy`, :math:`xz`, :math:`yz`).  The "thermo" data follows the
format of the *yaml* thermo style.  The "keywords" tag lists the
per-atom properties contained in the "data" columns, which contain a
list with one line per atom.  The keywords may be renamed using the
dump_modify command same as for the *custom* dump style.

.. code-block:: yaml

   ---
   creator: LAMMPS
   timestep: 0
   units: lj
   time: 0
   natoms: 4000
   boundary: [ p, p, p, p, p, p, ]
   thermo:
     - keywords: [ Step, Temp, E_pair, E_mol, TotEng, Press, ]
     - data: [ 0, 0, -27093.472213010766, 0, 0, 0, ]
   box:
     - [ 0, 16.795961913825074 ]
     - [ 0, 16.795961913825074 ]
     - [ 0, 16.795961913825074 ]
     - [ 0, 0, 0 ]
   keywords: [ id, type, x, y, z, vx, vy, vz, ix, iy, iz,  ]
   data:
     - [     1 , 1 ,  0.000000e+00 ,  0.000000e+00 ,  0.000000e+00 ,  -1.841579e-01 , -9.710036e-01 , -2.934617e+00 , 0 , 0 , 0, ]
     - [     2 , 1 ,  8.397981e-01 ,  8.397981e-01 ,  0.000000e+00 ,  -1.799591e+00 ,  2.127197e+00 ,  2.298572e+00 , 0 , 0 , 0, ]
     - [     3 , 1 ,  8.397981e-01 ,  0.000000e+00 ,  8.397981e-01 ,  -1.807682e+00 , -9.585130e-01 ,  1.605884e+00 , 0 , 0 , 0, ]

     [...]
   ...
   ---
   timestep: 100
   units: lj
   time: 0.5

     [...]

   ...

----------

Frequency of dump output:

Dumps are performed on timesteps that are a multiple of :math:`N`
(including timestep 0) and on the last timestep of a minimization if
the minimization converges.  Note that this means a dump will not be
performed on the initial timestep after the dump command is invoked,
if the current timestep is not a multiple of :math:`N`.  This behavior
can be changed via the :doc:`dump_modify first <dump_modify>` command,
which can also be useful if the dump command is invoked after a
minimization ended on an arbitrary timestep.

The value of :math:`N` can be changed between runs by using the
:doc:`dump_modify every <dump_modify>` command (not allowed for *dcd*
style).  The :doc:`dump_modify every <dump_modify>` command also
allows a variable to be used to determine the sequence of timesteps on
which dump files are written.  In this mode a dump on the first
timestep of a run will also not be written unless the
:doc:`dump_modify first <dump_modify>` command is used.

If you instead want to dump snapshots based on simulation time (in
time units of the :doc:`units` command), the :doc:`dump_modify
every/time <dump_modify>` command can be used.  This can be useful
when the timestep size varies during a simulation run, e.g. by use of
the :doc:`fix dt/reset <fix_dt_reset>` command.

----------

Dump filenames:

The specified dump filename determines how the dump file(s) is
written.  The default is to write one large text file, which is opened
when the dump command is invoked and closed when an :doc:`undump
<undump>` command is used or when LAMMPS exits.  For the *dcd* and
*xtc* styles, this is a single large binary file.

Many of the styles allow dump filenames to contain either or both of
two wildcard characters.  If a "\*" character appears in the filename,
then one file per snapshot is written and the "\*" character is
replaced with the timestep value.  For example, tmp.dump.\* becomes
tmp.dump.0, tmp.dump.10000, tmp.dump.20000, etc.  This option is not
available for the *dcd* and *xtc* styles.  Note that the
:doc:`dump_modify pad <dump_modify>` command can be used to ensure all
timestep numbers are the same length (e.g., 00010), which can make it
easier to read a series of dump files in order with some
post-processing tools.

If a "%" character appears in the filename, then each of P processors
writes a portion of the dump file, and the "%" character is replaced
with the processor ID from :math:`0` to :math:`P-1`.  For example,
tmp.dump.% becomes tmp.dump.0, tmp.dump.1, ... tmp.dump.:math:`P-1`,
etc.  This creates smaller files and can be a fast mode of output on
parallel machines that support parallel I/O for output. This option is
**not** available for the *dcd*, *xtc*, *xyz*, *grid/vtk*, and *yaml*
styles.

By default, :math:`P` is the the number of processors, meaning one file per
processor, but :math:`P` can be set to a smaller value via the *nfile* or
*fileper* keywords of the :doc:`dump_modify <dump_modify>` command.
These options can be the most efficient way of writing out dump files
when running on large numbers of processors.

Note that using the "\*" and "%" characters together can produce a
large number of small dump files!

.. deprecated:: 21Nov2023

The MPIIO package and the the corresponding "/mpiio" dump styles, except
for the unrelated "netcdf/mpiio" style were removed from LAMMPS.

----------

Compression of dump file data:

If the specified filename ends with ".bin" or ".lammpsbin", the dump
file (or files, if "\*" or "%" is also used) is written in binary
format.  A binary dump file will be about the same size as a text
version, but will typically write out much faster.  Of course, when
post-processing, you will need to convert it back to text format (see
the :ref:`binary2txt tool <binary>`) or write your own code to read
the binary file.  The format of the binary file can be understood by
looking at the :file:`tools/binary2txt.cpp` file.  This option is only
available for the *atom* and *custom* styles.

If the filename ends with ".gz", the dump file (or files, if "\*" or "%"
is also used) is written in gzipped format.  A gzipped dump file will be
about :math:`3\times` smaller than the text version, but will also take
longer to write.  This option is not available for the *dcd* and *xtc*
styles.

Note that styles that end with *gz* are identical in command syntax to
the corresponding styles without "gz", however, they generate
compressed files using the zlib library. Thus the filename suffix
".gz" is mandatory. This is an alternative approach to writing
compressed files via a pipe, as done by the regular dump styles, which
may be required on clusters where the interface to the high-speed
network disallows using the fork() library call (which is needed for a
pipe).  For the remainder of this page, you should thus consider the
*atom* and *atom/gz* styles (etc.) to be inter-changeable, with the
exception of the required filename suffix.

Similarly, styles that end with *zstd* are identical to the gz styles,
but use the Zstd compression library instead and require a ".zst"
suffix. See the :doc:`dump_modify <dump_modify>` page for details on
how to control the compression level in both variants.

----------

General triclinic simulation box output for the *atom* and *custom* styles:

As mentioned above, the simulation box can be defined as a general
triclinic box, which means that 3 arbitrary box edge vectors **A**,
**B**, **C** can be specified.  See the :doc:`Howto triclinic
<Howto_triclinic>` doc page for a detailed description of general
triclinic boxes.

This option is provided as a convenience for users who may be
converting data from solid-state crystallographic representations or
from DFT codes for input to LAMMPS.  However, as explained on the
:doc:`Howto_triclinic <Howto_triclinic>` doc page, internally, LAMMPS
only uses restricted triclinic simulation boxes.  This means the box
and per-atom information (e.g. coordinates, velocities) LAMMPS stores
are converted (rotated) from general to restricted triclinic form when
the system is created.

For dump output, if the :doc:`dump_modify triclinic/general
<dump_modify>` command is used, the box description and per-atom
coordinates and other per-atom vectors will be converted (rotated)
from restricted to general form when each dump file snapshots is
output.  This option can only be used if the simulation box was
initially created as general triclinic.  If the option is not used,
and the simulation box is general triclinic, then the dump file
snapshots will reflect the internal restricted triclinic geometry.

The dump_modify triclinic/general option affects 3 aspects of the dump
file output.

First, the format for the BOX BOUNDS is as follows

.. parsed-literal::

   ITEM: BOX BOUNDS abc origin
   ax ay az originx
   bx by bz originy
   cx cy cz originz

where the **A** edge vector of the box is (ax,ay,az) and similarly
for **B** and **C**.  The origin of all 3 edge vectors is (originx,
originy, originz).

Second, the coordinates of each atom are converted (rotated) so that
the atom is inside (or near) the general triclinic box defined by the
**A**, **B**, **C** edge vectors.  For style *atom*, this only alters
output for unscaled atom coords, via the :doc:`dump_modify scaled no
<dump_modify>` setting. For style *custom*, this alters output for
either unscaled or unwrapped output of atom coords, via the *x,y,z* or
*xu,yu,zu* attributes.  For output of scaled atom coords by both
styles, there is no difference between restricted and general
triclinic values.

Third, the output for any attribute of the *custom* style which
represents a per-atom vector quantity will be converted (rotated) to
be oriented consistent with the general triclinic box and its
orientation relative to the standard xyz coordinate axes.

This applies to the following *custom* style attributes:

* vx,vy,vz = atom velocities
* fx,fy,fz = forces on atoms
* mux,muy,muz = orientation of dipole moment of atom
* omegax,omegay,omegaz = angular velocity of spherical particle
* angmomx,angmomy,angmomz = angular momentum of aspherical particle
* tqx,tqy,tqz = torque on finite-size particles

For example, if the velocity of an atom in a restricted triclinic box
is along the x-axis, then it will be output for a general triclinic
box as a vector along the **A** edge vector of the box.

.. note::

   For style *custom*, the :doc:`dump_modify thresh <dump_modify>`
   command may access per-atom attributes either directly or
   indirectly through a compute or variable.  If the attribute is an
   atom coordinate or one of the vectors mentioned above, its value
   will *NOT* be a general triclinic (rotated) value.  Rather it will
   be a restricted triclinic value.

----------

Arguments for different styles:

The sections below describe per-atom, local, and per grid cell
attributes which can be used as arguments to the various styles.

Note that in the discussion below, for styles which can reference
values from a compute or fix or custom atom property, like the
*custom*\ , *cfg*\ , *grid* or *local* styles, the bracketed index
:math:`i` can be specified using a wildcard asterisk with the index to
effectively specify multiple values.  This takes the form "\*" or
"\*n" or "m\*" or "m\*n".  If :math:`N` is the number of columns in
the array, then an asterisk with no numeric values means all column
indices from 1 to :math:`N`.  A leading asterisk means all indices
from 1 to n (inclusive).  A trailing asterisk means all indices from m
to :math:`N` (inclusive).  A middle asterisk means all indices from m
to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  For example, these two dump commands are
equivalent, since the :doc:`compute stress/atom <compute_stress_atom>`
command creates a per-atom array with six columns:

.. code-block:: LAMMPS

   compute myPress all stress/atom NULL
   dump 2 all custom 100 tmp.dump id myPress[*]
   dump 2 all custom 100 tmp.dump id myPress[1] myPress[2] myPress[3] &
                                     myPress[4] myPress[5] myPress[6]

----------

Per-atom attributes used as arguments to the *custom* and *cfg* styles:

The *id*, *mol*, *proc*, *procp1*, *type*, *element*, *mass*, *vx*,
*vy*, *vz*, *fx*, *fy*, *fz*, *q* attributes are self-explanatory.

*Id* is the atom ID.  *Mol* is the molecule ID, included in the data
file for molecular systems.  *Proc* is the ID of the processor (0 to
:math:`N_\text{procs}-1`) that currently owns the atom.  *Procp1* is the
proc ID+1, which can be convenient in place of a *type* attribute (1 to
:math:`N_\text{types}`) for coloring atoms in a visualization program.
*Type* is the atom type (1 to :math:`N_\text{types}`).  *Element* is
typically the chemical name of an element, which you must assign to each
type via the :doc:`dump_modify element <dump_modify>` command.  More
generally, it can be any string you wish to associated with an atom
type.  *Mass* is the atom mass. The quantities *vx*, *vy*, *vz*, *fx*,
*fy*, *fz*, and *q* are components of atom velocity and force and atomic
charge.

There are several options for outputting atom coordinates.  The *x*,
*y*, and *z* attributes write atom coordinates "unscaled", in the
appropriate distance :doc:`units <units>` (:math:`\AA`,
:math:`\sigma`, etc.).  Use *xs*, *ys*, and *zs* if you want the
coordinates "scaled" to the box size so that each value is 0.0 to 1.0.
If the simulation box is triclinic (tilted), then all atom coords will
still be between 0.0 and 1.0.  The actual unscaled :math:`(x,y,z)`
coordinate is :math:`x_s a + y_s b + z_s c`, where :math:`(a,b,c)` are
the non-orthogonal vectors of the simulation box edges, as discussed on
the :doc:`Howto triclinic <Howto_triclinic>` page.

Use *xu*, *yu*, and *zu* if you want the coordinates "unwrapped" by the
image flags for each atom.  Unwrapped means that if the atom has passed
through a periodic boundary one or more times, the value is printed for
what the coordinate would be if it had not been wrapped back into the
periodic box.  Note that using *xu*, *yu*, and *zu* means that the
coordinate values may be far outside the box bounds printed with the
snapshot.  Using *xsu*, *ysu*, and *zsu* is similar to using *xu*, *yu*,
and *zu*, except that the unwrapped coordinates are scaled by the box
size. Atoms that have passed through a periodic boundary will have the
corresponding coordinate increased or decreased by 1.0.

The image flags can be printed directly using the *ix*, *iy*, and *iz*
attributes.  For periodic dimensions, they specify which image of the
simulation box the atom is considered to be in.  An image of 0 means
it is inside the box as defined.  A value of 2 means add 2 box lengths
to get the true value.  A value of :math:`-1` means subtract 1 box length to
get the true value.  LAMMPS updates these flags as atoms cross
periodic boundaries during the simulation.

The *mux*, *muy*, and *muz* attributes are specific to dipolar systems
defined with an atom style of *dipole*\ .  They give the orientation of
the atom's point dipole moment.  The *mu* attribute gives the magnitude
of the atom's dipole moment.

The *radius* and *diameter* attributes are specific to spherical
particles that have a finite size, such as those defined with an atom
style of *sphere*\ .

The *omegax*, *omegay*, and *omegaz* attributes are specific to
finite-size spherical particles that have an angular velocity.  Only
certain atom styles, such as *sphere*, define this quantity.

The *angmomx*, *angmomy*, and *angmomz* attributes are specific to
finite-size aspherical particles that have an angular momentum.  Only
the *ellipsoid* atom style defines this quantity.

The *tqx*, *tqy*, and *tqz* attributes are for finite-size particles
that can sustain a rotational torque due to interactions with other
particles.

The *c_ID* and *c_ID[I]* attributes allow per-atom vectors or arrays
calculated by a :doc:`compute <compute>` to be output.  The ID in the
attribute should be replaced by the actual ID of the compute that has
been defined previously in the input script.  See the :doc:`compute
<compute>` command for details.  There are computes for calculating the
per-atom energy, stress, centro-symmetry parameter, and coordination
number of individual atoms.

Note that computes which calculate global or local quantities, as
opposed to per-atom quantities, cannot be output in a dump custom
command.  Instead, global quantities can be output by the
:doc:`thermo_style custom <thermo_style>` command, and local quantities
can be output by the dump local command.

If *c_ID* is used as a attribute, then the per-atom vector calculated by
the compute is printed.  If *c_ID[i]* is used, then :math:`i` must be in
the range from 1 to :math:`M`, which will print the :math:`i`\ th column
of the per-atom array with :math:`M` columns calculated by the compute.
See the discussion above for how :math:`i` can be specified with a
wildcard asterisk to effectively specify multiple values.

The *f_ID* and *f_ID[I]* attributes allow vector or array per-atom
quantities calculated by a :doc:`fix <fix>` to be output.  The ID in the
attribute should be replaced by the actual ID of the fix that has been
defined previously in the input script.  The :doc:`fix ave/atom
<fix_ave_atom>` command is one that calculates per-atom quantities.
Since it can time-average per-atom quantities produced by any
:doc:`compute <compute>`, :doc:`fix <fix>`, or atom-style :doc:`variable
<variable>`, this allows those time-averaged results to be written to a
dump file.

If *f_ID* is used as a attribute, then the per-atom vector calculated by
the fix is printed.  If *f_ID[i]* is used, then :math:`i` must be in the
range from 1 to :math:`M`, which will print the :math:`i`\ th column of
the per-atom array with :math:`M` columns calculated by the fix.  See
the discussion above for how :math:`i` can be specified with a wildcard
asterisk to effectively specify multiple values.

The *v_name* attribute allows per-atom vectors calculated by a
:doc:`variable <variable>` to be output.  The name in the attribute
should be replaced by the actual name of the variable that has been
defined previously in the input script.  Only an atom-style variable can
be referenced, since it is the only style that generates per-atom
values.  Variables of style *atom* can reference individual atom
attributes, per-atom attributes, thermodynamic keywords, or invoke other
computes, fixes, or variables when they are evaluated, so this is a very
general means of creating quantities to output to a dump file.

The *i_name*, *d_name*, *i2_name*, *d2_name* attributes refer to
custom per-atom integer and floating-point vectors or arrays that have
been added via the :doc:`fix property/atom <fix_property_atom>`
command.  When that command is used specific names are given to each
attribute which are the "name" portion of these keywords.  For arrays
*i2_name* and *d2_name*, the column of the array must also be included
following the name in brackets (e.g., d2_xyz[i], i2_mySpin[i], where
:math:`i` is in the range from 1 to :math:`M`, where :math:`M` is the
number of columns in the custom array).  See the discussion above for
how :math:`i` can be specified with a wildcard asterisk to effectively
specify multiple values.

See the :doc:`Modify <Modify>` page for information on how to add
new compute and fix styles to LAMMPS to calculate per-atom quantities
which could then be output into dump files.

----------

Attributes used as arguments to the *local* style:

The *index* attribute can be used to generate an index number from 1
to N for each line written into the dump file, where N is the total
number of local datums from all processors, or lines of output that
will appear in the snapshot.  Note that because data from different
processors depend on what atoms they currently own, and atoms migrate
between processor, there is no guarantee that the same index will be
used for the same info (e.g. a particular bond) in successive
snapshots.

The *c_ID* and *c_ID[I]* attributes allow local vectors or arrays
calculated by a :doc:`compute <compute>` to be output.  The ID in the
attribute should be replaced by the actual ID of the compute that has
been defined previously in the input script.  See the
:doc:`compute <compute>` command for details.  There are computes for
calculating local information such as indices, types, and energies for
bonds and angles.

Note that computes which calculate global or per-atom quantities, as
opposed to local quantities, cannot be output in a dump local command.
Instead, global quantities can be output by the :doc:`thermo_style
custom <thermo_style>` command, and per-atom quantities can be output
by the dump custom command.

If *c_ID* is used as a attribute, then the local vector calculated by
the compute is printed.  If *c_ID[I]* is used, then I must be in the
range from 1-M, which will print the Ith column of the local array
with M columns calculated by the compute.  See the discussion above
for how I can be specified with a wildcard asterisk to effectively
specify multiple values.

The *f_ID* and *f_ID[I]* attributes allow local vectors or arrays
calculated by a :doc:`fix <fix>` to be output.  The ID in the attribute
should be replaced by the actual ID of the fix that has been defined
previously in the input script.

If *f_ID* is used as a attribute, then the local vector calculated by
the fix is printed.  If *f_ID[I]* is used, then I must be in the
range from 1-M, which will print the Ith column of the local with M
columns calculated by the fix.  See the discussion above for how I can
be specified with a wildcard asterisk to effectively specify multiple
values.

Here is an example of how to dump bond info for a system, including
the distance and energy of each bond:

.. code-block:: LAMMPS

   compute 1 all property/local batom1 batom2 btype
   compute 2 all bond/local dist eng
   dump 1 all local 1000 tmp.dump index c_1[1] c_1[2] c_1[3] c_2[1] c_2[2]

----------

Attributes used as arguments to the *grid* and *grid/vtk* styles:

The attributes that begin with *c_ID* and *f_ID* both take
colon-separated fields *gname* and *dname*.  These refer to a grid
name and data field name which is defined by the compute or fix.  Note
that a compute or fix can define one or more grids (of different
sizes) and one or more data fields for each of those grids.  The sizes
of all grids output in a single dump grid command must be the same.

The *c_ID:gname:dname* and *c_ID:gname:dname[I]* attributes allow
per-grid vectors or arrays calculated by a :doc:`compute <compute>` to
be output.  The ID in the attribute should be replaced by the actual
ID of the compute that has been defined previously in the input
script.

If *c_ID:gname:dname* is used as a attribute, then the per-grid vector
calculated by the compute is printed.  If *c_ID:gname:dname[I]* is
used, then I must be in the range from 1-M, which will print the Ith
column of the per-grid array with M columns calculated by the compute.
See the discussion above for how I can be specified with a wildcard
asterisk to effectively specify multiple values.

The *f_ID:gname:dname* and *f_ID:gname:dname[I]* attributes allow
per-grid vectors or arrays calculated by a :doc:`fix <fix>` to be
output.  The ID in the attribute should be replaced by the actual ID
of the fix that has been defined previously in the input script.

If *f_ID:gname:dname* is used as a attribute, then the per-grid vector
calculated by the fix is printed.  If *f_ID:gname:dname[I]* is used,
then I must be in the range from 1-M, which will print the Ith column
of the per-grid with M columns calculated by the fix.  See the
discussion above for how I can be specified with a wildcard asterisk
to effectively specify multiple values.

----------

Restrictions
""""""""""""

To write gzipped dump files, you must either compile LAMMPS with the
-DLAMMPS_GZIP option or use the styles from the COMPRESS package.
See the :doc:`Build settings <Build_settings>` page for details.

While a dump command is active (i.e., has not been stopped by using
the :doc:`undump command <undump>`), no commands may be used that will
change the timestep (e.g., :doc:`reset_timestep <reset_timestep>`).
LAMMPS will terminate with an error otherwise.

The *atom/gz*, *cfg/gz*, *custom/gz*, and *xyz/gz* styles are part of
the COMPRESS package.  They are only enabled if LAMMPS was built with
that package.  See the :doc:`Build package <Build_package>` page for
more info.

The *xtc*, *dcd*, and *yaml* styles are part of the EXTRA-DUMP package.
They are only enabled if LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`dump atom/adios <dump_adios>`, :doc:`dump custom/adios <dump_adios>`,
:doc:`dump cfg/uef <dump_cfg_uef>`, :doc:`dump h5md <dump_h5md>`,
:doc:`dump image <dump_image>`, :doc:`dump molfile <dump_molfile>`,
:doc:`dump netcdf <dump_netcdf>`, :doc:`dump netcdf/mpiio <dump_netcdf>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`,
:doc:`write_dump <write_dump>`

Default
"""""""

The defaults for the *image* and *movie* styles are listed on the
:doc:`dump image <dump_image>` page.
