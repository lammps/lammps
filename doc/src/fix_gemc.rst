.. index:: fix gemc

fix gemc command
================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID gemc N M X V T displace maxdlogvolratio seed keyword ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* gemc = style name of this fix command
* N = invoke this fix every N steps
* M = average number of translations (and rotations of molecules) to attempt every N steps
* X = average number of atom or molecule exchanges between the two boxes to attempt every N steps
* V = average number of volume exchanges between the two boxes to attempt every N steps
* T = temperature of the Gibbs ensemble (temperature units)
* displace = maximum Monte Carlo translation distance (distance units)
* maxdlogvolratio = maximum change of ln(V1/V2) in a volume exchange (unitless)
* seed = random # seed (positive integer)
* zero or more keywords may be appended
* keyword = *full_energy* or *tune* or *mol* or *maxangle*

  .. parsed-literal::

       *full_energy* = compute the full energy of the system for every move
       *mol* value = template-ID
         template-ID = ID of molecule template specified in a separate :doc:`molecule <molecule>` command
       *maxangle* value = maximum rotation angle of molecules (degrees)
       *tune* values = Nt Atrans Avol
         Nt = adjust step sizes every Nt invocations of this fix
         Atrans = target acceptance ratio of translations and rotations (0 < Atrans < 1)
         Avol = target acceptance ratio of volume exchanges (0 < Avol < 1)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all gemc 1 100 20 2 0.9 0.3 0.05 29494
   fix mc all gemc 10 1000 200 10 120.0 1.0 0.1 4711 full_energy
   fix mc all gemc 1 100 100 2 0.9 0.1 0.1 29494 tune 20 0.4 0.4
   fix mc all gemc 1 100 100 2 250.0 0.5 0.1 4711 mol co2mol maxangle 30

Description
"""""""""""

.. versionadded:: 4Jul2026

This fix performs Monte Carlo (MC) moves in the Gibbs ensemble (GEMC)
using two separate simulation boxes that are in thermodynamic contact
without having an interface between them :ref:`(Panagiotopoulos)
<Panagiotopoulos1>`, :ref:`(Frenkel) <Frenkel3>`.  The boxes exchange
volume, so that the pressure is the same in both, and they exchange
atoms, so that the chemical potential is the same in both, while the
total volume, the total number of atoms, and the temperature *T* are
fixed.  The most common application is the calculation of vapor-liquid
coexistence: when started at a total density inside the two-phase
region, one box ends up containing the vapor and the other box the
liquid, and the average densities of the boxes are the two coexisting
densities at temperature *T*.  The :doc:`Gibbs ensemble Howto
<Howto_gemc>` explains how to set up and analyze such a calculation
step by step.

The two boxes are two :doc:`partitions <Run_options>` of LAMMPS, so
LAMMPS must be started with exactly two partitions, for example
``mpirun -np 2 lmp -in in.gemc -partition 1 1``.  Both partitions read
the same input script.  Commands that should differ between the boxes,
like the initial box size or number of atoms, can use
:doc:`world-style variables <variable>` or the :doc:`partition
<partition>` command.  The two partitions can use different numbers of
processors, which is useful because the liquid box usually contains
many more atoms than the vapor box, e.g. ``-partition 1 3`` on 4
processors.  The fix command itself must be the same in both
partitions, except for the *displace* and *maxangle* values and the
*full_energy* keyword, which may differ between the boxes.  The fix
group, the names of all groups, the number of atom types, the atom
style, and the molecule template must also be the same, since atoms
are exchanged between the boxes; otherwise the fix stops with an
error at the start of a run.

Every *N* timesteps the fix performs a total of *M* + *X* + *V* MC
moves.  The type of each move is chosen randomly with probabilities
proportional to *M*, *X*, and *V*, so these are the average numbers of
the three kinds of moves.  Both boxes always attempt the same kind of
move at the same time:

* A *translation* moves a randomly chosen atom of the fix group in each
  box by a random displacement inside a sphere of radius *displace*.
  The radius is limited to half of the smallest box width, and when a
  box runs on more than one processor, to half of the smallest subdomain
  width, since an atom can only move to a neighboring subdomain during
  the move.
  The two boxes accept or reject their translations independently with
  the Metropolis criterion :math:`\min[1, \exp(-\Delta U/k_B T)]`.

* An *exchange* removes a randomly chosen atom of the fix group from one
  box, the donor box, and inserts it at a random position in the other
  box.  The donor box is chosen randomly with equal probability.  The
  inserted atom keeps the atom type, charge, group membership, and
  velocity of the removed atom, so that time integration continues
  consistently in hybrid MD/MC simulations.  The move is accepted with probability

  .. math::

     \min\left[1, \frac{N_d\,V_r}{(N_r+1)\,V_d}
     \exp\left(-\frac{\Delta U_d + \Delta U_r}{k_B T}\right)\right]

  where :math:`N_d, V_d` and :math:`N_r, V_r` are the number of atoms
  in the fix group and the volume of the donor and receiving box before
  the move.

* A *volume exchange* changes :math:`\ln(V_1/V_2)` by a random amount
  between -*maxdlogvolratio* and +*maxdlogvolratio* while keeping the
  total volume :math:`V_1+V_2` constant.  Each box is scaled
  uniformly with its lower corner kept fixed, and all atom coordinates
  are scaled with it (for molecules, see below).  The move is accepted with probability

  .. math::

     \min\left[1, \left(\frac{V_1'}{V_1}\right)^{N_1+1}
     \left(\frac{V_2'}{V_2}\right)^{N_2+1}
     \exp\left(-\frac{\Delta U_1 + \Delta U_2}{k_B T}\right)\right]

  where :math:`N_1, N_2` are the total numbers of atoms in the boxes.

Exchange and volume moves are always accepted or rejected by both boxes
together.  The energy *U* is the total potential energy of a box, as
computed by the *thermo_pe* compute, see below.

.. versionchanged:: TBD

An exchange keeps the atom type, charge, group membership, and velocity
of the removed atom; previously an atom of type 1 with a random velocity
was inserted.  Atom styles with per-atom masses now stop with an error
instead of a warning.

Atoms that are not in the fix group are never translated or
exchanged, but they interact with the other atoms and are scaled by
volume moves.  Different atom types can be used, e.g. to study the
coexistence of a mixture; the exchange move then transfers the type
and charge of the randomly chosen atom.

.. versionadded:: TBD

**Molecules:** With the *mol* keyword, the fix moves and exchanges whole
molecules instead of atoms.  All atoms of the fix group must then belong
to molecules that are copies of the molecule template given by
*template-ID*: each molecule must have the same number of atoms with the
same atom types as the template, and consecutive atom IDs in the order
of the template atoms.  This is the case for molecules created with the
:doc:`create_atoms <create_atoms>` command using the same template, and
for all molecules inserted by this fix.  Each molecule must also have
the same number of bonds, angles, dihedrals, and impropers as the
template.  Without the *mol* keyword, the atoms of the fix group must not
be part of bonds, angles, dihedrals, or impropers.  The template provides the
bond topology of inserted molecules; the atom style must allow it, and
enough space for bonds, angles, and special neighbors must be reserved,
e.g. with the *extra/special/per/atom* keyword of :doc:`create_box
<create_box>` or :doc:`read_data <read_data>`.  With molecules:

* Each translation move is, with equal probability, a translation of a
  randomly chosen molecule by a random displacement inside a sphere of
  radius *displace*, or a rotation of a randomly chosen molecule about
  its center of mass around a random axis by a random angle between
  -*maxangle* and +*maxangle*.  The default *maxangle* is 30 degrees.
  Rotations that would move an atom farther than the limit for the
  translation radius are rejected.
  These are rigid-body moves; to sample the internal degrees of
  freedom of flexible molecules, combine the fix with time integration
  (see below).

* An exchange removes a randomly chosen molecule from the donor box and
  inserts it with a random orientation at a random position in the
  other box.  The inserted molecule keeps the conformation, the charges,
  the group membership, and the velocities (rotated with the molecule)
  of its atoms, so the move is also correct for flexible molecules.
  The acceptance probability is the one given above, with :math:`N`
  the number of molecules.  For the reason given below, an exchange
  into a box narrower than twice the size of the molecule is rejected.

* A volume exchange moves the center of mass of each molecule with the
  box, while the shape of the molecules is not changed.  :math:`N` in
  the acceptance probability is the number of molecules plus the number
  of atoms that do not belong to a molecule.  Since bonded interactions
  use the closest periodic image of bond partners, volume exchanges
  that would make a box narrower than twice the size of its largest
  molecule are rejected.  This applies to all systems with molecule IDs,
  also without the *mol* keyword.

Moves of molecules always use the total energy of the system, see
below, including the bonded interactions.

.. versionadded:: TBD

**Charges:** Charged atoms and molecules, and long-range solvers set
with :doc:`kspace_style <kspace_style>`, are supported.  With a
long-range solver, the total energy is computed for all moves.  An
exchange moves the charge of an atom or of the atoms of a molecule from
one box to the other box.  If the exchanged atoms or molecules carry a
net charge, the boxes are not charge neutral, and the fix prints a
warning.  For physically meaningful results, exchange only neutral
molecules, or neutral atoms.

.. versionadded:: TBD

**Triclinic boxes:** The boxes may be triclinic.  A volume exchange
scales the box lengths and tilt factors by the same factor, so the box
shape does not change.

Choosing the parameters: *displace* is usually chosen so that roughly
30% to 50% of the translations in the liquid box are accepted, and
*maxdlogvolratio* so that roughly 30% to 50% of the volume exchanges
are accepted.  The *tune* keyword can adjust both automatically, see
below.  The acceptance ratio of exchanges is usually low for dense
liquids (a few percent or less), so *X* must be large enough to
exchange every atom several times during the run.  The acceptance
counts are available as output of this fix, see below.

.. versionadded:: TBD

The *tune* keyword adjusts *displace*, *maxdlogvolratio*, and, for
molecules, *maxangle* during the run.  Every *Nt* invocations of the
fix, the acceptance ratio of the translations, rotations, and volume
exchanges attempted since the last adjustment is compared to the
targets *Atrans* and *Avol*, and the step size is multiplied by the
ratio of the measured to the target acceptance ratio, limited to the
range 0.5 to 1.5.  The maximum displacement and rotation angle are
adjusted separately for each box, since the vapor box accepts much
larger moves than the liquid box.  The displacement is limited to half
the smallest box width, to half the smallest subdomain width when a
box runs on more than one processor, and, when single-atom energies are
used (see below), to half the distance by which the ghost atoms extend
beyond the pair cutoff, so that the ghost atoms rarely need to be
rebuilt.  The maximum change of
ln(V1/V2) is limited to 1.0, and the maximum rotation angle to 180
degrees.  An adjustment is only made after at least 20 moves of the
respective kind were attempted.

.. note::

   Changing the step sizes based on the history of the simulation
   violates detailed balance, so the *tune* keyword should only be used
   during equilibration.  For the production run, re-define the fix
   without the *tune* keyword using the adjusted step sizes, which are
   available as elements 7, 8, and 11 of the output vector.  Since the
   maximum displacement differs between the boxes, it can be passed
   through a variable that is evaluated in each partition:

   .. code-block:: LAMMPS

      fix             mc all gemc 1 100 100 2 0.9 0.1 0.1 29494 tune 20 0.4 0.4
      run             10000
      variable        disp equal $(f_mc[7])
      variable        dlv equal $(f_mc[8])
      unfix           mc
      fix             mc all gemc 1 100 100 2 0.9 ${disp} ${dlv} 29494
      run             100000

If the fix is used together with time integration, e.g. :doc:`fix nvt
<fix_nh>`, a hybrid MD/MC simulation is performed.  In this case the
thermostat temperature should be the same as *T*, and the temperature
compute used by the thermostat should account for the changing number
of atoms, for example:

.. code-block:: LAMMPS

   compute mdtemp all temp
   compute_modify mdtemp dynamic/dof yes
   fix mdnvt all nvt temp 0.9 0.9 0.5
   fix_modify mdnvt temp mdtemp

In pure MC simulations without time integration, the atom velocities
are irrelevant, and LAMMPS will print a warning that atoms will not
move, which can be ignored.  In both cases, use :doc:`compute_modify
thermo_temp dynamic/dof yes <compute_modify>` so that the thermodynamic
output uses the current number of atoms in each box.

.. versionadded:: TBD

By default, translation and exchange moves compute only the change of
the pair energy of the moved atom, which is much faster than computing
the total energy of the box: the cost of such a move does not depend on
the number of atoms.  This is only possible when the pair style supports
the *single()* function, is not a many-body potential, does not use
:doc:`pair_modify tail yes <pair_modify>`, and when no fix contributes
to the potential energy, no molecules are exchanged (*mol* keyword), the
exchanged atoms are not charged, and no long-range solver is used.
*displace* must also not be larger than the distance by which the
:doc:`ghost atoms <Developer_par_comm>` extend beyond the pair cutoff,
which is the neighbor skin distance set by the :doc:`neighbor
<neighbor>` command, or larger with :doc:`comm_modify cutoff
<comm_modify>`.  If any of these conditions is not met, a warning is
printed and the total energy is computed for every move.  The same is
done for a move in a box whose width is smaller than the pair cutoff.
Volume moves always compute the total energy.  The *full_energy*
keyword requests that the total energy is computed for all moves.

.. versionchanged:: TBD

The cost of translations and exchanges with single-atom energies no
longer grows with the number of atoms in the box.  Ghost atoms are
rebuilt only when an atom would move farther than the distance above
outside of the subdomain of its processor, or beyond the neighboring
subdomains, so a larger communication cutoff makes this less frequent.

Some fixes have an associated potential energy.  Examples of such fixes
include: :doc:`efield <fix_efield>`, :doc:`gravity <fix_gravity>`,
:doc:`addforce <fix_addforce>`, :doc:`restrain <fix_restrain>`, and
:doc:`wall fixes <fix_wall>`.  For that energy to be included in the
total potential energy of the system (the quantity used by the MC
moves), you MUST enable the :doc:`fix_modify <fix_modify>` *energy*
option for that fix.

Neighbor lists are rebuilt every *N* timesteps that this fix is
invoked, so you should not set *N* too small when combining the fix
with time integration.  In pure MC simulations *N* = 1 with a large
number of moves per invocation is most efficient.

During a run, the fix prints a progress report to the screen and to
the log file of the whole run (not to those of the partitions) every 1%
of the run.  It lists the numbers of accepted and attempted moves since
the last report (translations of the first box, exchanges and volume
changes of both boxes) and the volume, number of atoms, and number
density of each box.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the fix to :doc:`binary restart files
<restart>`.  This includes information about the random number
generator seeds, the next timestep for MC moves, and the numbers of MC
move attempts and successes.  Each partition writes its own restart
file, so the restart file name must be different for the two
partitions, e.g. by using a world-style variable.  See the
:doc:`read_restart <read_restart>` command for info on how to
re-specify a fix in an input script that reads a restart file, so that
the operation of the fix continues in an uninterrupted fashion.

.. note::

   For this to work correctly, the timestep must **not** be changed
   after reading the restart with :doc:`reset_timestep
   <reset_timestep>`.  The fix will try to detect it and stop with an
   error.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

.. versionchanged:: TBD

This fix computes a global vector of length 11, which can be accessed by
various :doc:`output commands <Howto_output>`, e.g. as *f_ID[1]* in
the :doc:`thermo_style <thermo_style>` command.  The vector values are
the following cumulative counts and current settings for the box of the
partition where they are accessed:

  #. translation attempts
  #. translation successes
  #. exchange attempts
  #. exchange successes
  #. volume change attempts
  #. volume change successes
  #. current maximum translation distance *displace* of this box
  #. current maximum change of ln(V1/V2) *maxdlogvolratio*
  #. rotation attempts (molecules only)
  #. rotation successes (molecules only)
  #. current maximum rotation angle *maxangle* of this box (degrees)

The translation counts include only translations, not rotations.

The vector values calculated by this fix are "intensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

This fix requires exactly two partitions, a 3d simulation with a box
that is periodic in all three dimensions, and atom IDs.  Atom styles
with per-atom masses or with per-atom properties other than charge,
molecule ID, and bond topology (e.g. point dipoles, orientations, or
custom properties from :doc:`fix property/atom <fix_property_atom>`)
are not supported, since exchanged atoms would lose these properties.

With the *mol* keyword, only molecules of a single kind (one molecule
template) can be moved and exchanged.  Atom style *template* is not
supported.  Constraints with :doc:`fix rigid <fix_rigid>`, :doc:`fix
shake <fix_shake>`, or :doc:`fix rattle <fix_shake>` are not supported;
molecules stay rigid in pure MC simulations, since all their moves are
rigid-body moves.

Inserted atoms and molecules get new atom and molecule IDs that are
larger than all IDs in their box.  The IDs of removed atoms and
molecules are therefore not reused, except when the one with the largest
ID was removed, and the largest ID keeps growing during a long
simulation.  This increases the memory used by the atom map, and
the IDs may eventually exceed the maximum allowed value.  Use
:doc:`atom_modify map hash <atom_modify>` to limit the memory use.
When atoms (not molecules) are exchanged, the :doc:`reset_atoms id
<reset_atoms>` command can be used between runs to compress the IDs.

Do not set :doc:`neigh_modify once yes <neigh_modify>` or else this fix
will never be called.  Reneighboring is **required**.

Use of multiple *fix gemc* commands in the same input script can be
problematic.

Related commands
""""""""""""""""

:doc:`fix gcmc <fix_gcmc>`,
:doc:`fix widom <fix_widom>`,
:doc:`fix atom/swap <fix_atom_swap>`,
:doc:`partition <partition>`

Defaults
""""""""

By default, single-atom energies are used when possible, i.e. the
*full_energy* keyword is not set, and the step sizes are not adjusted,
i.e. the *tune* keyword is not set.  Atoms are exchanged (no *mol*
keyword), and maxangle = 30 degrees.

----------

.. _Panagiotopoulos1:

**(Panagiotopoulos)** Panagiotopoulos, Mol Phys, 61, 813-826 (1987).

.. _Frenkel3:

**(Frenkel)** Frenkel and Smit, Understanding Molecular Simulation,
3rd edition, Academic Press, London, 2023.
