.. index:: fix tmd

fix tmd command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID tmd rho_final file1 N file2

* ID, group-ID are documented in :doc:`fix <fix>` command
* tmd = style name of this fix command
* rho_final = desired value of rho at the end of the run (distance units)
* file1 = filename to read target structure from
* N = dump TMD statistics every this many timesteps, 0 = no dump
* file2 = filename to write TMD statistics to (only needed if N > 0)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve
   fix 2 tmdatoms tmd 1.0 target_file 100 tmd_dump_file

Description
"""""""""""

Perform targeted molecular dynamics (TMD) on a group of atoms.  A
holonomic constraint is used to force the atoms to move towards (or
away from) the target configuration.  The parameter "rho" is
monotonically decreased (or increased) from its initial value to
rho_final at the end of the run.

Rho has distance units and is a measure of the root-mean-squared
distance (RMSD) between the current configuration of the atoms in the
group and the target coordinates listed in file1.  Thus a value of
rho_final = 0.0 means move the atoms all the way to the final
structure during the course of the run.

The target file1 can be ASCII text or a gzipped text file (detected by
a .gz suffix).  The format of the target file1 is as follows:

.. parsed-literal::

   0.0 25.0 xlo xhi
   0.0 25.0 ylo yhi
   0.0 25.0 zlo zhi
   125     24.97311   1.69005     23.46956 0 0 -1
   126     1.94691    2.79640     1.92799  1 0 0
   127     0.15906    3.46099     0.79121  1 0 0
   ...

The first 3 lines may or may not be needed, depending on the format of
the atoms to follow.  If image flags are included with the atoms, the
1st 3 lo/hi lines must appear in the file.  If image flags are not
included, the 1st 3 lines should not appear.  The 3 lines contain the
simulation box dimensions for the atom coordinates, in the same format
as in a LAMMPS data file (see the :doc:`read_data <read_data>` command).

The remaining lines each contain an atom ID and its target x,y,z
coordinates.  The atom lines (all or none of them) can optionally be
followed by 3 integer values: nx,ny,nz.  For periodic dimensions, they
specify which image of the box the atom is considered to be in, i.e. a
value of N (positive or negative) means add N times the box length to
the coordinate to get the true value.

The atom lines can be listed in any order, but every atom in the group
must be listed in the file.  Atoms not in the fix group may also be
listed; they will be ignored.

TMD statistics are written to file2 every N timesteps, unless N is
specified as 0, which means no statistics.

The atoms in the fix tmd group should be integrated (via a fix nve,
nvt, npt) along with other atoms in the system.

Restarts can be used with a fix tmd command.  For example, imagine a
10000 timestep run with a rho_initial = 11 and a rho_final = 1.  If a
restart file was written after 2000 time steps, then the configuration
in the file would have a rho value of 9.  A new 8000 time step run
could be performed with the same rho_final = 1 to complete the
conformational change at the same transition rate.  Note that for
restarted runs, the name of the TMD statistics file should be changed
to prevent it being overwritten.

For more information about TMD, see :ref:`(Schlitter1) <Schlitter1>` and
:ref:`(Schlitter2) <Schlitter2>`.

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.

This fix can ramp its rho parameter over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

All TMD fixes must be listed in the input script after all integrator
fixes (nve, nvt, npt) are applied.  This ensures that atoms are moved
before their positions are corrected to comply with the constraint.

Atoms that have a TMD fix applied should not be part of a group to
which a SHAKE fix is applied.  This is because LAMMPS assumes there
are not multiple competing holonomic constraints applied to the same
atoms.

To read gzipped target files, you must compile LAMMPS with the
-DLAMMPS_GZIP option.  See the :doc:`Build settings <Build_settings>`
doc page for details.

**Related commands:** none

**Default:** none

----------

.. _Schlitter1:

**(Schlitter1)** Schlitter, Swegat, Mulders, "Distance-type reaction
coordinates for modelling activated processes", J Molecular Modeling,
7, 171-177 (2001).

.. _Schlitter2:

**(Schlitter2)** Schlitter and Klahn, "The free energy of a reaction
coordinate at multiple constraints: a concise formulation", Molecular
Physics, 101, 3439-3443 (2003).
