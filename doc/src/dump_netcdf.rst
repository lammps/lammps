.. index:: dump netcdf

dump netcdf command
===================

dump netcdf/mpiio command
=========================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID netcdf N file args
   dump ID group-ID netcdf/mpiio N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be imaged
* *netcdf* or *netcdf/mpiio*  = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = list of atom attributes, same as for :doc:`dump_style custom <dump>`

Examples
""""""""

.. parsed-literal::

   dump 1 all netcdf 100 traj.nc type x y z vx vy vz
   dump_modify 1 append yes at -1 thermo yes
   dump 1 all netcdf/mpiio 1000 traj.nc id type x y z
   dump 1 all netcdf 1000 traj.\*.nc id type x y z

Description
"""""""""""

Dump a snapshot of atom coordinates every N timesteps in Amber-style
NetCDF file format.  NetCDF files are binary, portable and
self-describing.  This dump style will write only one file on the root
node.  The dump style *netcdf* uses the `standard NetCDF library <netcdf-home_>`_.  All data is collected on one processor and then
written to the dump file.  Dump style *netcdf/mpiio* uses the
`parallel NetCDF library <pnetcdf-home_>`_ and MPI-IO to write to the dump
file in parallel; it has better performance on a larger number of
processors.  Note that style *netcdf* outputs all atoms sorted by atom
tag while style *netcdf/mpiio* outputs atoms in order of their MPI
rank.

NetCDF files can be directly visualized via the following tools:

Ovito (http://www.ovito.org/). Ovito supports the AMBER convention and
all extensions of this dump style.

* VMD (http://www.ks.uiuc.edu/Research/vmd/).
* AtomEye (http://www.libatoms.org/). The libAtoms version of AtomEye
  contains a NetCDF reader that is not present in the standard
  distribution of AtomEye.

In addition to per-atom data, :doc:`thermo <thermo>` data can be included in the
dump file. The data included in the dump file is identical to the data specified
by :doc:`thermo_style <thermo_style>`.

.. _netcdf-home: http://www.unidata.ucar.edu/software/netcdf/

.. _pnetcdf-home: http://trac.mcs.anl.gov/projects/parallel-netcdf/

----------

Restrictions
""""""""""""

The *netcdf* and *netcdf/mpiio* dump styles are part of the
USER-NETCDF package.  They are only enabled if LAMMPS was built with
that package. See the :doc:`Build package <Build_package>` doc page for
more info.

----------

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
