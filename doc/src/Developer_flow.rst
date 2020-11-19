How a timestep works
--------------------

The first and most fundamental operation within LAMMPS to understand is
how a timestep is structured.  Timestepping is performed by calling
methods of the Integrate class instance within the Update class.  Since
Integrate is a base class, it will point to an instance of a derived
class corresponding to what is selected by the :doc:`run_style
<run_style>` input script command.

In this section, the timestep implemented by the Verlet class is
described.  A similar timestep protocol is implemented by the Respa
class, for the r-RESPA hierarchical timestepping method.

The Min base class performs energy minimization, so does not perform a
literal timestep.  But it has logic similar to what is described here,
to compute forces and invoke fixes at each iteration of a minimization.
Differences between time integration and minimization are highlighted at
the end of this section.

The Verlet class is encoded in the ``src/verlet.cpp`` and ``verlet.h``
files.  It implements the velocity-Verlet timestepping algorithm.  The
workhorse method is ``Verlet::run()``, but first we highlight several
other methods in the class.

- The ``init()`` method is called at the beginning of each dynamics
  run.  It simply sets some internal flags, based on user settings in
  other parts of the code.

- The ``setup()`` or ``setup_minimal()`` methods are also called before
  each run.  The velocity-Verlet method requires current forces be
  calculated before the first timestep, so these routines compute
  forces due to all atomic interactions, using the same logic that
  appears in the timestepping described next.  A few fixes are also
  invoked, using the mechanism described in the next section.  Various
  counters are also initialized before the run begins.  The
  ``setup_minimal()`` method is a variant that has a flag for performing
  less setup.  This is used when runs are continued and information
  from the previous run is still valid.  For example, if repeated
  short LAMMPS runs are being invoked, interleaved by other commands,
  via the *pre no* and *every* options of the run command, the
  ``setup_minimal()`` method is used.

- The ``force_clear()`` method initializes force and other arrays to
  zero before each timestep, so that forces (torques, etc) can be
  accumulated.

Now for the ``Verlet::run()`` method.  Its basic structure in hi-level pseudo
code is shown below.  In the actual code in ``src/verlet.cpp`` some of
these operations are conditionally invoked.

.. code-block:: python

   loop over N timesteps:
     if timeout condition: break
     ev_set()

     fix->initial_integrate()
     fix->post_integrate()

     nflag = neighbor->decide()
     if nflag:
       fix->pre_exchange()
       domain->pbc()
       domain->reset_box()
       comm->setup()
       neighbor->setup_bins()
       comm->exchange()
       comm->borders()
       fix->pre_neighbor()
       neighbor->build()
       fix->post_neighbor()
     else:
       comm->forward_comm()

     force_clear()
     fix->pre_force()

     pair->compute()
     bond->compute()
     angle->compute()
     dihedral->compute()
     improper->compute()
     kspace->compute()

     fix->pre_reverse()
     comm->reverse_comm()

     fix->post_force()
     fix->final_integrate()
     fix->end_of_step()

     if any output on this step:
       output->write()

   # after loop
   fix->post_run()


The ``ev_set()`` method (in the parent Integrate class), sets two flags
(*eflag* and *vflag*) for energy and virial computation.  Each flag
encodes whether global and/or per-atom energy and virial should be
calculated on this timestep, because some fix or variable or output will
need it.  These flags are passed to the various methods that compute
particle interactions, so that they either compute and tally the
corresponding data or can skip the extra calculations if the energy and
virial are not needed.  See the comments for the ``Integrate::ev_set()``
method which document the flag values.

At various points of the timestep, fixes are invoked,
e.g. ``fix->initial_integrate()``.  In the code, this is actually done
via the Modify class which stores all the Fix objects and lists of which
should be invoked at what point in the timestep.  Fixes are the LAMMPS
mechanism for tailoring the operations of a timestep for a particular
simulation.  As described elsewhere, each fix has one or more methods,
each of which is invoked at a specific stage of the timestep, as show in
the timestep pseudo-code.  All the active fixes defined in an input
script, that are flagged to have an ``initial_integrate()`` method are
invoked at the beginning of each timestep.  Examples are :doc:`fix nve
<fix_nve>` or :doc:`fix nvt or fix npt <fix_nh>` which perform the
start-of-timestep velocity-Verlet integration operations to update
velocities by a half-step, and coordinates by a full step.  The
``post_integrate()`` method is next for operations that need to happen
immediately after those updates.  Only a few fixes use this, e.g. to
reflect particles off box boundaries in the :doc:`FixWallReflect class
<fix_wall_reflect>`.

The ``decide()`` method in the Neighbor class determines whether
neighbor lists need to be rebuilt on the current timestep (conditions
can be changed using the :doc:`neigh_modify every/delay/check
<neigh_modify>` command.  If not, coordinates of ghost atoms are
acquired by each processor via the ``forward_comm()`` method of the Comm
class.  If neighbor lists need to be built, several operations within
the inner if clause of the pseudo-code are first invoked.  The
``pre_exchange()`` method of any defined fixes is invoked first.
Typically this inserts or deletes particles from the system.

Periodic boundary conditions are then applied by the Domain class via
its ``pbc()`` method to remap particles that have moved outside the
simulation box back into the box.  Note that this is not done every
timestep, but only when neighbor lists are rebuilt.  This is so that
each processor's sub-domain will have consistent (nearby) atom
coordinates for its owned and ghost atoms.  It is also why dumped atom
coordinates may be slightly outside the simulation box if not dumped
on a step where the neighbor lists are rebuilt.

The box boundaries are then reset (if needed) via the ``reset_box()``
method of the Domain class, e.g. if box boundaries are shrink-wrapped to
current particle coordinates.  A change in the box size or shape
requires internal information for communicating ghost atoms (Comm class)
and neighbor list bins (Neighbor class) be updated.  The ``setup()``
method of the Comm class and ``setup_bins()`` method of the Neighbor
class perform the update.

The code is now ready to migrate atoms that have left a processor's
geometric sub-domain to new processors.  The ``exchange()`` method of
the Comm class performs this operation.  The ``borders()`` method of the
Comm class then identifies ghost atoms surrounding each processor's
sub-domain and communicates ghost atom information to neighboring
processors.  It does this by looping over all the atoms owned by a
processor to make lists of those to send to each neighbor processor.  On
subsequent timesteps, the lists are used by the ``Comm::forward_comm()``
method.

Fixes with a ``pre_neighbor()`` method are then called.  These typically
re-build some data structure stored by the fix that depends on the
current atoms owned by each processor.

Now that each processor has a current list of its owned and ghost
atoms, LAMMPS is ready to rebuild neighbor lists via the ``build()``
method of the Neighbor class.  This is typically done by binning all
owned and ghost atoms, and scanning a stencil of bins around each
owned atom's bin to make a Verlet list of neighboring atoms within the
force cutoff plus neighbor skin distance.

In the next portion of the timestep, all interaction forces between
particles are computed, after zeroing the per-atom force vector via the
``force_clear()`` method.  If the newton flag is set to *on* by the
newton command, forces are added to both owned and ghost atoms, otherwise
only to owned (aka local) atoms.

Pairwise forces are calculated first, which enables the global virial
(if requested) to be calculated cheaply (at O(N) cost instead of O(N**2)
at the end of the ``Pair::compute()`` method), by a dot product of atom
coordinates and forces.  By including owned and ghost atoms in the dot
product, the effect of periodic boundary conditions is correctly
accounted for.  Molecular topology interactions (bonds, angles,
dihedrals, impropers) are calculated next (if supported by the current
atom style).  The final contribution is from long-range Coulombic
interactions, invoked by the KSpace class.

The ``pre_reverse()`` method in fixes is used for operations that have to
be done *before* the upcoming reverse communication (e.g. to perform
additional data transfers or reductions for data computed during the
force computation and stored with ghost atoms).

If the newton flag is on, forces on ghost atoms are communicated and
summed back to their corresponding owned atoms.  The ``reverse_comm()``
method of the Comm class performs this operation, which is essentially
the inverse operation of sending copies of owned atom coordinates to
other processor's ghost atoms.

At this point in the timestep, the total force on each (local) atom is
known.  Additional force constraints (external forces, SHAKE, etc) are
applied by Fixes that have a ``post_force()`` method.  The second half
of the velocity-Verlet integration, ``final_integrate()`` is then
performed (another half-step update of the velocities) via fixes like
nve, nvt, npt.

At the end of the timestep, fixes that contain an ``end_of_step()``
method are invoked.  These typically perform a diagnostic calculation,
e.g. the ave/time and ave/spatial fixes.  The final operation of the
timestep is to perform any requested output, via the ``write()`` method
of the Output class.  There are 3 kinds of LAMMPS output: thermodynamic
output to the screen and log file, snapshots of atom data to a dump
file, and restart files.  See the :doc:`thermo_style <thermo_style>`,
:doc:`dump <dump>`, and :doc:`restart <restart>` commands for more
details.

The the flow of control during energy minimization iterations is
similar to that of a molecular dynamics timestep.  Forces are computed,
neighbor lists are built as needed, atoms migrate to new processors, and
atom coordinates and forces are communicated to neighboring processors.
The only difference is what Fix class operations are invoked when.  Only
a subset of LAMMPS fixes are useful during energy minimization, as
explained in their individual doc pages.  The relevant Fix class methods
are ``min_pre_exchange()``, ``min_pre_force()``, and ``min_post_force()``.
Each fix is invoked at the appropriate place within the minimization
iteration.  For example, the ``min_post_force()`` method is analogous to
the ``post_force()`` method for dynamics; it is used to alter or constrain
forces on each atom, which affects the minimization procedure.

After all iterations are completed there is a ``cleanup`` step which
calls the ``post_run()`` method of fixes to perform operations only required
at the end of a calculations (like freeing temporary storage or creating
final outputs).
