.. index:: fix gemc

fix gemc command
================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID gemc N M X V T displace maxvol seed

* ID, group-ID are documented in :doc:`fix <fix>` command
* gemc = style name of this fix command
* N = invoke this fix every N steps
* M = average number of MC moves to attempt every N steps
* X = average number of GEMC exchanges to attempt every N steps
* V = average number of GEMC volume changes to attempt every N steps
* T = temperature of the Gibbs ensemble (temperature units)
* displace = maximum Monte Carlo translation distance (length units)
* maxdlogvolratio = maximum change in ln(vol1/vol2) (unitless)
* seed = random # seed (positive integer)
* zero keyword/value pairs may be appended to args

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all gemc 100 10 20 30 1.1 0.5 100.0 29494

Description
"""""""""""

.. versionadded:: 4Jul2026

This fix performs Gibbs ensemble Monte Carlo (GEMC) exchanges of atoms
and volume between two simulation cells at specified *T*.  It also
attempts Monte Carlo (MC) moves (atom translations) within the
simulation cell.  This is usually used to establish thermodynamic
equilibrium between bulk vapor and liquid phases, as discussed in
:ref:`(Frenkel) <Frenkel3>`.  If used with the :doc:`fix nvt <fix_nh>`
command, hybrid MD/MC simulations in the Gibbs ensemble (equal pressure,
equal chemical potential, constant total volume, and constant
temperature) can be performed.  Specific uses include computing
vapor-liquid coexistence curves.

Every *N* timesteps the fix attempts GEMC atom exchanges, GEMC volume
changes, and MC moves of atoms.  On those timesteps, the average number
of attempted GEMC atom exchanges is *X*, the average number of volume
changes is *V*, and the average number of attempted MC moves is *M*.

This fix requires that LAMMPS be run with two partitions that
instantiate the two simulation cells.  This requires using the
:doc:`-partition command-line switch <Run_options>`.  For example, on a
workstation with 12 cores, ``-partition 2x6`` could be used.  For better
performance, it is recommended that the two partitions be initialized at
two different densities.  This allows assigning only one core to the
partition running the vapor phase and all the other cores to the
partition running the liquid phase.  This can be achieved in a single
script by using the :doc:`partition <partition>` command to initiate the
two partitions at different densities. On a workstation with 12
processor cores, the command line option ``-partition 1 8`` can be used
to assign 1 core to the first (vapor) partition and 8 cores to the
second (liquid) partition.

If used with :doc:`fix nvt <fix_nh>`, the temperature of the Gibbs
ensemble, *T*, should be set to be equivalent to the target temperature
used in fix nvt.  Otherwise, the imaginary reservoir will not be in
thermal equilibrium with the simulation cell.  Also, it is important
that the temperature used by *fix nvt* is dynamically updated, which can
be achieved as follows:

.. code-block:: LAMMPS

   compute mdtemp mdatoms temp
   compute_modify mdtemp dynamic/dof yes
   fix mdnvt mdatoms nvt temp 300.0 300.0 10.0
   fix_modify mdnvt temp mdtemp

Note that neighbor lists are re-built every timestep that this fix is
invoked, so you should not set *N* to be too small.  However, periodic
rebuilds are necessary in order to avoid dangerous rebuilds and missed
interactions.  Specifically, avoid performing so many MC translations
per timestep that atoms can move beyond the neighbor list skin distance.
See the :doc:`neighbor <neighbor>` command for details.

When an atom is inserted in either partition, its coordinates are chosen
at a random position within the current simulation cell, and new atom
velocities are randomly chosen from the specified temperature
distribution given by *T*.

Some fixes have an associated potential energy. Examples of such fixes
include: :doc:`efield <fix_efield>`, :doc:`gravity <fix_gravity>`,
:doc:`addforce <fix_addforce>`, :doc:`langevin <fix_langevin>`,
:doc:`restrain <fix_restrain>`, :doc:`temp/berendsen
<fix_temp_berendsen>`, :doc:`temp/rescale <fix_temp_rescale>`, and
:doc:`wall fixes <fix_wall>`.  For that energy to be included in the
total potential energy of the system (the quantity used when performing
GEMC atom exchange, GEMC volume exchange and MC moves), you MUST enable
the :doc:`fix_modify <fix_modify>` *energy* option for that fix.  The
doc pages for individual :doc:`fix <fix>` commands specify if this
should be done.

Use of this fix typically will cause the number of atoms in each cell to
fluctuate, therefore, you will want to use the :doc:`compute_modify
dynamic/dof <compute_modify>` command to ensure that the current number
of atoms is used as a normalizing factor each time temperature is
computed. A simple example of this is:

.. code-block:: LAMMPS

   compute_modify thermo_temp dynamic/dof yes

A more complicated example is listed earlier on this page in the context
of NVT dynamics.

.. note::

   If the density of the cell is initially very small or zero, and
   increases to a much larger density after a period of equilibration,
   then certain quantities that are only calculated once at the start
   (kspace parameters) may no longer be accurate.  The solution is to
   start a new simulation after the equilibrium density has been
   reached.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the fix to :doc:`binary restart files
<restart>`.  This includes information about the random number generator
seed, the next timestep for MC exchanges, the number of MC step attempts
and successes etc.  See the :doc:`read_restart <read_restart>` command
for info on how to re-specify a fix in an input script that reads a
restart file, so that the operation of the fix continues in an
uninterrupted fashion.

.. note::

   For this to work correctly, the timestep must **not** be changed
   after reading the restart with :doc:`reset_timestep
   <reset_timestep>`.  The fix will try to detect it and stop with an
   error.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix computes a global vector of length 8, which can be accessed by
various :doc:`output commands <Howto_output>`.  The vector values are
the following global cumulative quantities:

  #. translation attempts
  #. translation successes
  #. exchange attempts
  #. exchange successes
  #. volume change attempts
  #. volume change successes

The vector values calculated by this fix are "intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Do not set :doc:`neigh_modify once yes <neigh_modify>` or else this fix
will never be called.  Reneighboring is **required**.

*Fix gemc* currently **only** supports MC moves and exchanges on
individual atoms.

Use of multiple *fix gemc* commands in the same input script can be
problematic.

Related commands
""""""""""""""""

:doc:`fix gcmc <fix_gcmc>`,
:doc:`fix widom <fix_widom>`,
:doc:`fix atom/swap <fix_atom_swap>`

Defaults
""""""""

None

----------

.. _Frenkel3:

**(Frenkel)** Frenkel and Smit, Understanding Molecular Simulation,
Academic Press, London, 2002.
