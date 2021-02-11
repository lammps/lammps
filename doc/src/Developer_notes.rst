Notes for developers and code maintainers
-----------------------------------------

This section documents how some of the code functionality within
LAMMPS works at a conceptual level.  Comments on code in source files
typically document what a variable stores, what a small section of
code does, or what a function does and its input/outputs.  The topics
on this page are intended to document code functionality at a higher level.

Fix contributions to instantaneous energy, virial, and cumulative energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixes can calculate contributions to the instantaneous energy and/or
virial of the system, both in a global and peratom sense.  Fixes that
perform thermostatting or barostatting can calculate the cumulative
energy they add to or subtract from the system, which is accessed by
the *ecouple* and *econserve* thermodynamic keywords.  This subsection
explains how both work and what flags to set in a new fix to enable
this functionality.

Let's start with thermostatting and barostatting fixes.  Examples are
the :doc:`fix langevin <fix_langevin>` and :doc:`fix npt <fix_nh>`
commands.  Here is what the fix needs to do:

* Set the variable *ecouple_flag* = 1 in the constructor.  Also set
  *scalar_flag* = 1, *extscalar* = 1, and *global_freq* to a timestep
  increment which matches how often the fix is invoked.
* Implement a compute_scalar() method that returns the cumulative
  energy added or subtracted by the fix, e.g. by rescaling the
  velocity of atoms.  The sign convention is that subtracted energy is
  positive, added energy is negative.  This must be the total energy
  added to the entire system, i.e. an "extensive" quantity, not a
  per-atom energy.  Cumulative means the summed energy since the fix
  was instantiated, even across multiple runs.  This is because the
  energy is used by the *econserve* thermodynamic keyword to check
  that the fix is conserving the total energy of the system,
  i.e. potential energy + kinetic energy + coupling energy = a
  constant.

And here is how the code operates:

* The Modify class makes a list of all fixes that set *ecouple_flag* = 1.
* The :doc:`thermo_style custom <thermo_style>` command defines
  *ecouple* and *econserve* keywords.
* These keywords sum the energy contributions from all the
  *ecouple_flag* = 1 fixes by invoking the energy_couple() method in
  the Modify class, which calls the compute_scalar() method of each
  fix in the list.

------------------

Next, here is how a fix contributes to the instantaneous energy and
virial of the system.  First, it sets any or all of these flags to a
value of 1 in their constructor:

* *energy_global_flag* to contribute to global energy, example: :doc:`fix indent <fix_indent>`
* *energy_peratom_flag* to contribute to peratom energy, :doc:`fix cmap <fix_cmap>`
* *virial_global_flag* to contribute to global virial, example: :doc:`fix wall <fix_wall>`
* *virial_peratom_flag* to contribute to peratom virial, example: :doc:`fix wall <fix_wall>`

The fix must also do the following:

* For global energy, implement a compute_scalar() method that returns
  the energy added or subtracted on this timestep.  Here the sign
  convention is that added energy is positive, subtracted energy is
  negative.
* For peratom energy, invoke the ev_init(eflag,vflag) function each
  time the fix is invoked, which initializes per-atom energy storage.
  The value of eflag may need to be stored from an earlier call to the
  fix during the same timestep.  See how the :doc:`fix cmap
  <fix_cmap>` command does this in src/MOLECULE/fix_cmap.cpp.  When an
  energy for one or more atoms is calculated, invoke the ev_tally()
  function to tally the contribution to each atom.  Both the ev_init()
  and ev_tally() methods are in the parent Fix class.
* For global and/or peratom virial, invoke the v_init(vflag) function
  each time the fix is invoked, which initializes virial storage.
  When forces on one or more atoms are calculated, invoke the
  v_tally() function to tally the contribution.  Both the v_init() and
  v_tally() methods are in the parent Fix class.  Note that there are
  several variants of v_tally(); choose the one appropriate to your
  fix.

.. note::

   The ev_init() and ev_tally() methods also account for global and
   peratom virial contributions.  Thus you do not need to invoke the
   v_init() and v_tally() methods, if the fix also calculates peratom
   energies.

The fix must also specify whether (by default) to include or exclude
these contributions to the global/peratom energy/virial of the system.
For the fix to include the contributions, set either of both of these
variables in the constructor:

* *thermo_energy* = 1, for global and peratom energy
* *thermo_virial* = 1, for global and peratom virial

Note that these variables are zeroed in fix.cpp.  Thus if you don't
set the variables, the contributions will be excluded (by default)

However, the user has ultimate control over whether to include or
exclude the contributions of the fix via the :doc:`fix modify
<fix_modify>` command:

* fix modify *energy yes* to include global and peratom energy contributions
* fix modify *virial yes* to include global and peratom virial contributions

If the fix contributes to any of the global/peratom energy/virial
values for the system, it should be explained on the fix doc page,
along with the default values for the *energy yes/no* and *virial
yes/no* settings of the :doc:`fix modify <fix_modify>` command.

Finally, these 4 contributions are included in the output of 4
computes:

* global energy in :doc:`compute pe <compute_pe>`
* peratom energy in :doc:`compute pe/atom <compute_pe_atom>`
* global virial in :doc:`compute pressure <compute_pressure>`
* peratom virial in :doc:`compute stress/atom <compute_stress_atom>`

These computes invoke a method of the Modify class to include
contributions from fixes that have the corresponding flags set,
e.g. *energy_peratom_flag* and *thermo_energy* for :doc:`compute
pe/atom <compute_pe_atom>`.

Note that each compute has an optional keyword to either include or
exclude all contributions from fixes.  Also note that :doc:`compute pe
<compute_pe>` and :doc:`compute pressure <compute_pressure>` are what
is used (by default) by :doc:`thermodynamic output <thermo_style>` to
calculate values for its *pe* and *press* keywords.

KSpace PPPM FFT grids
^^^^^^^^^^^^^^^^^^^^^

The various :doc:`KSpace PPPM <kspace_style>` styles in LAMMPS use
FFTs to solve Poisson's equation.  This subsection describes:

* how FFT grids are defined
* how they are decomposed across processors
* how they are indexed by each processor
* how particle charge and electric field values are mapped to/from
  the grid

An FFT grid cell is a 3d volume; grid points are corners of a grid
cell and the code stores values assigned to grid points in vectors or
3d arrays.  A global 3d FFT grid has points indexed 0 to N-1 inclusive
in each dimension.

Each processor owns two subsets of the grid, each subset is
brick-shaped.  Depending on how it is used, these subsets are
allocated as a 1d vector or 3d array.  Either way, the ordering of
values within contiguous memory x fastest, then y, z slowest.

For the ``3d decomposition`` of the grid, the global grid is
partitioned into bricks that correspond to the sub-domains of the
simulation box that each processor owns.  Often, this is a regular 3d
array (Px by Py by Pz) of bricks, where P = number of processors =
Px * Py * Pz.  More generally it can be a tiled decomposition, where
each processor owns a brick and the union of all the bricks is the
global grid.  Tiled decompositions are produced by load balancing with
the RCB algorithm; see the :doc:`balance rcb <balance>` command.

For the ``FFT decompostion`` of the grid, each processor owns a brick
that spans the entire x dimension of the grid while the y and z
dimensions are partitioned as a regular 2d array (P1 by P2), where P =
P1 * P2.

The following indices store the inclusive bounds of the brick a
processor owns, within the global grid:

.. parsed-literal::

   nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in = 3d decomposition brick
   nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft = FFT decomposition brick
   nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out = 3d decomposition brick + ghost cells

The ``in`` and ``fft`` indices are from 0 to N-1 inclusive in each
dimension, where N is the grid size.

The ``out`` indices index an array which stores the ``in`` subset of
the grid plus ghost cells that surround it.  These indices can thus be
< 0 or >= N.

The number of ghost cells a processor owns in each of the 6 directions
is a function of:

.. parsed-literal::

   neighbor skin distance (since atoms can move outside a proc subdomain)
   qdist = offset or charge from atom due to TIP4P fictitious charge
   order = mapping stencil size
   shift = factor used when order is an even number (see below)

Here is an explanation of how the PPPM variables ``order``,
``nlower`` / ``nupper``, ``shift``, and ``OFFSET`` work. They are the
relevant variables that determine how atom charge is mapped to grid
points and how field values are mapped from grid points to atoms:

.. parsed-literal::

   order = # of nearby grid points in each dim that atom charge/field are mapped to/from
   nlower,nupper = extent of stencil around the grid point an atom is assigned to
   OFFSET = large integer added/subtracted when mapping to avoid int(-0.75) = 0 when -1 is the desired result

The particle_map() method assigns each atom to a grid point.

If order is even, say 4:

.. parsed-literal::

   atom is assigned to grid point to its left (in each dim)
   shift = OFFSET
   nlower = -1, nupper = 2, which are offsets from assigned grid point
   window of mapping grid pts is thus 2 grid points to left of atom, 2 to right

If order is odd, say 5:

.. parsed-literal::

   atom is assigned to left/right grid pt it is closest to (in each dim)
   shift = OFFSET + 0.5
   nlower = 2, nupper = 2
   if point is in left half of cell, then window of affected grid pts is 3 grid points to left of atom, 2 to right
   if point is in right half of cell, then window of affected grid pts is 2 grid points to left of atom, 3 to right

These settings apply to each dimension, so that if order = 5, an
atom's charge is mapped to 125 grid points that surround the atom.
