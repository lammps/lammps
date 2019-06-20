.. index:: rerun

rerun command
=============

Syntax
""""""


.. parsed-literal::

   rerun file1 file2 ... keyword args ...

* file1,file2,... = dump file(s) to read
* one or more keywords may be appended, keyword *dump* must appear and be last
  
  .. parsed-literal::
  
     keyword = *first* or *last* or *every* or *skip* or *start* or *stop* or *dump*
      *first* args = Nfirst
        Nfirst = dump timestep to start on
      *last* args = Nlast
        Nlast = dumptimestep to stop on
      *every* args = Nevery
        Nevery = read snapshots matching every this many timesteps
      *skip* args = Nskip
        Nskip = read one out of every Nskip snapshots
      *start* args = Nstart
        Nstart = timestep on which pseudo run will start
      *stop* args = Nstop
        Nstop = timestep to which pseudo run will end
      *dump* args = same as :doc:`read_dump <read_dump>` command starting with its field arguments



Examples
""""""""


.. parsed-literal::

   rerun dump.file dump x y z vx vy vz
   rerun dump1.txt dump2.txt first 10000 every 1000 dump x y z
   rerun dump.vels dump x y z vx vy vz box yes format molfile lammpstrj
   rerun dump.dcd dump x y z box no format molfile dcd
   rerun ../run7/dump.file.gz skip 2 dump x y z box yes

Description
"""""""""""

Perform a pseudo simulation run where atom information is read one
snapshot at a time from a dump file(s), and energies and forces are
computed on the shapshot to produce thermodynamic or other output.

This can be useful in the following kinds of scenarios, after an
initial simulation produced the dump file:

* Compute the energy and forces of snaphots using a different potential.
* Calculate one or more diagnostic quantities on the snapshots that
  weren't computed in the initial run.  These can also be computed with
  settings not used in the initial run, e.g. computing an RDF via the
  :doc:`compute rdf <compute_rdf>` command with a longer cutoff than was
  used initially.
* Calculate the portion of per-atom forces resulting from a subset of
  the potential.  E.g. compute only Coulombic forces.  This can be done
  by only defining only a Coulombic pair style in the rerun script.
  Doing this in the original script would result in different (bad)
  dynamics.


Conceptually, using the rerun command is like running an input script
that has a loop in it (see the :doc:`next <next>` and :doc:`jump <jump>`
commands).  Each iteration of the loop reads one snapshot from the
dump file via the :doc:`read\_dump <read_dump>` command, sets the
timestep to the appropriate value, and then invokes a :doc:`run <run>`
command for zero timesteps to simply compute energy and forces, and
any other :doc:`thermodynamic output <thermo_style>` or diagnostic info
you have defined.  This computation also invokes any fixes you have
defined that apply constraints to the system, such as :doc:`fix shake <fix_shake>` or :doc:`fix indent <fix_indent>`.

Note that a simulation box must already be defined before using the
rerun command.  This can be done by the :doc:`create\_box <create_box>`,
:doc:`read\_data <read_data>`, or :doc:`read\_restart <read_restart>`
commands.

Also note that reading per-atom information from dump snapshots is
limited to the atom coordinates, velocities and image flags as
explained in the :doc:`read\_dump <read_dump>` command.  Other atom
properties, which may be necessary to compute energies and forces,
such as atom charge, or bond topology information for a molecular
system, are not read from (or even contained in) dump files.  Thus
this auxiliary information should be defined in the usual way, e.g. in
a data file read in by a :doc:`read\_data <read_data>` command, before
using the rerun command.

Also note that the frequency of thermodynamic or dump output from the
rerun simulation will depend on settings made in the rerun script, the
same as for output from any LAMMPS simulation.  See further info below
as to what that means if the timesteps for snapshots read from dump
files do not match the specified output frequency.


----------


If more than one dump file is specified, the dump files are read one
after the other.  It is assumed that snapshot timesteps will be in
ascending order.  If a snapshot is encountered that is not in
ascending order, it will skip the snapshot until it reads one that is.
This allows skipping of a duplicate snapshot (same timestep),
e.g. that appeared at the end of one file and beginning of the next.
However if you specify a series of dump files in an incorrect order
(with respect to the timesteps they contain), you may skip large
numbers of snapshots

Note that the dump files specified as part of the *dump* keyword can
be parallel files, i.e. written as multiple files either per processor
and/or per snapshot.  If that is the case they will also be read in
parallel which can make the rerun command operate dramatically faster
for large systems.  See the doc page for the `read\_dump <read_dump>`_ and
:doc:`dump <dump>` commands which describe how to read and write
parallel dump files.

The *first*\ , *last*\ , *every*\ , *skip* keywords determine which
snapshots are read from the dump file(s).  Snapshots are skipped until
they have a timestamp >= *Nfirst*\ .  When a snapshot with a timestamp >
*Nlast* is encountered, the rerun command finishes.  Note below that
the defaults for *first* and *last* are to read all snapshots.  If the
*every* keyword is set to a value > 0, then only snapshots with
timestamps that are a multiple of *Nevery* are read (the first
snapshot is always read).  If *Nevery* = 0, then this criterion is
ignored, i.e. every snapshot is read that meets the other criteria.
If the *skip* keyword is used, then after the first snapshot is read,
every Nth snapshot is read, where N = *Nskip*\ .  E.g. if *Nskip* = 3,
then only 1 out of every 3 snapshots is read, assuming the snapshot
timestamp is also consistent with the other criteria.

The *start* and *stop* keywords do not affect which snapshots are read
from the dump file(s).  Rather, they have the same meaning that they
do for the :doc:`run <run>` command.  They only need to be defined if
(a) you are using a :doc:`fix <fix>` command that changes some value
over time, and (b) you want the reference point for elapsed time (from
start to stop) to be different than the *first* and *last* settings.
See the doc page for individual fixes to see which ones can be used
with the *start/stop* keywords.  Note that if you define neither of
the *start*\ /\ *stop* or *first*\ /\ *last* keywords, then LAMMPS treats the
pseudo run as going from 0 to a huge value (effectively infinity).
This means that any quantity that a fix scales as a fraction of
elapsed time in the run, will essentially remain at its initial value.
Also note that an error will occur if you read a snapshot from the
dump file with a timestep value larger than the *stop* setting you
have specified.

The *dump* keyword is required and must be the last keyword specified.
Its arguments are passed internally to the :doc:`read\_dump <read_dump>`
command.  The first argument following the *dump* keyword should be
the *field1* argument of the :doc:`read\_dump <read_dump>` command.  See
the :doc:`read\_dump <read_dump>` doc page for details on the various
options it allows for extracting information from the dump file
snapshots, and for using that information to alter the LAMMPS
simulation.


----------


In general, a LAMMPS input script that uses a rerun command can
include and perform all the usual operations of an input script that
uses the :doc:`run <run>` command.  There are a few exceptions and
points to consider, as discussed here.

Fixes that perform time integration, such as :doc:`fix nve <fix_nve>` or
:doc:`fix npt <fix_nh>` are not invoked, since no time integration is
performed.  Fixes that perturb or constrain the forces on atoms will
be invoked, just as they would during a normal run.  Examples are :doc:`fix indent <fix_indent>` and :doc:`fix langevin <fix_langevin>`.  So you
should think carefully as to whether that makes sense for the manner
in which you are reprocessing the dump snapshots.

If you only want the rerun script to perform an analysis that does
not involve pair interactions, such as use compute msd to calculated
displacements over time, you do not need to define a :doc:`pair style <pair_style>`, which may also mean neighbor lists will not
need to be calculated which saves time.  The :doc:`comm\_modify cutoff <comm_modify>` command can also be used to insure ghost
atoms are acquired from far enough away for operations like bond and
angle evaluations, if no pair style is being used.

Every time a snapshot is read, the timestep for the simulation is
reset, as if the :doc:`reset\_timestep <reset_timestep>` command were
used.  This command has some restrictions as to what fixes can be
defined.  See its doc page for details.  For example, the :doc:`fix deposit <fix_deposit>` and :doc:`fix dt/reset <fix_dt_reset>` fixes
are in this category.  They also make no sense to use with a rerun
command.

If time-averaging fixes like :doc:`fix ave/time <fix_ave_time>` are
used, they are invoked on timesteps that are a function of their
*Nevery*\ , *Nrepeat*\ , and *Nfreq* settings.  As an example, see the
:doc:`fix ave/time <fix_ave_time>` doc page for details.  You must
insure those settings are consistent with the snapshot timestamps that
are read from the dump file(s).  If an averaging fix is not invoked on
a timestep it expects to be, LAMMPS will flag an error.

The various forms of LAMMPS output, as defined by the
:doc:`thermo\_style <thermo_style>`, :doc:`thermo <thermo>`,
:doc:`dump <dump>`, and :doc:`restart <restart>` commands occur with
specified frequency, e.g. every N steps.  If the timestep for a dump
snapshot is not a multiple of N, then it will be read and processed,
but no output will be produced.  If you want output for every dump
snapshot, you can simply use N=1 for an output frequency, e.g. for
thermodynamic output or new dump file output.


----------


Restrictions
""""""""""""


To read gzipped dump files, you must compile LAMMPS with the
-DLAMMPS\_GZIP option.  See the :doc:`Build settings <Build_settings>`
doc page for details.

Related commands
""""""""""""""""

:doc:`read\_dump <read_dump>`

Default
"""""""

The option defaults are first = 0, last = a huge value (effectively
infinity), start = same as first, stop = same as last, every = 0, skip
= 1;


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
