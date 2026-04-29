.. index:: fix msevb

fix msevb command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID msevb keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* msevb = style name of this fix command
* one or more keyword/value pairs must be specified

  .. parsed-literal::

     *coupling* values = style params
       style = *none* or *raiteri2011* or *vuilleumier1998* or *grimme2015*
         *none* = no off-diagonal coupling (diagonal EVB only)
         *raiteri2011* params = *lambda* V *zeta* Z
           V = coupling prefactor (energy units)
           Z = Gaussian decay parameter (inverse distance-squared units)
         *vuilleumier1998* params = *v12* V *alpha* A *gamma* G
           V = coupling prefactor (energy units)
           A = linear decay parameter (inverse distance units)
           G = Gaussian decay parameter (inverse distance-squared units)
         *grimme2015* params = *a* A *b* B
           A = coupling prefactor (energy units)
           B = energy decay parameter (inverse energy-squared units)
     *reaction* values = pre-ID post-ID mapfile cutoff [per-reaction keywords]
       pre-ID = molecule template ID for the pre-reaction topology
       post-ID = molecule template ID for the post-reaction topology
       mapfile = path to a REACTER-format map file
       cutoff = distance between initiator atoms at which reaction is considered (distance units)
       zero or more per-reaction keywords may follow in any order:
         *coupling* style params = per-reaction coupling override (same syntax as global *coupling*)
         *taper* T = distance at which coupling begins to taper to zero (distance units)
         *offset* E = energy offset added to this state's diagonal Hamiltonian element (energy units)
         *shells* N = maximum shell depth for this reaction (default: inherits global *shells*)
     *taper* T = global default taper start distance (distance units)
     *shells* value = N
       N = global default for maximum shell depth (default: 1)
     *fermi_dirac* value = T
       T = electronic temperature for Fermi-Dirac occupancy smearing (temperature units)
     *product_states* = enumerate product states for multi-site reactions (no value)
     *scf_topology* value = *yes* or *no*
       re-evaluate EVB after each permanent transfer until convergence (default: no)
     *scf_max_iter* value = N
       N = maximum SCF iterations per timestep (default: 10)
     *output* values = filename every
       filename = path to per-timestep EVB output file
       every = write output every this many steps

Examples
""""""""

.. code-block:: LAMMPS

   molecule pre  pre_h3o.mol
   molecule post post_h3o.mol
   fix evb all msevb &
     reaction pre post react.map 2.5 taper 2.0 shells 3 coupling grimme2015 a 0.6 b 0.25 &
     output msevb.out 1
   fix_modify evb energy yes

.. code-block:: LAMMPS

   molecule pre1 pre1.mol
   molecule post1 post1.mol
   molecule pre2 pre2.mol
   molecule post2 post2.mol
   fix evb all msevb &
     coupling raiteri2011 lambda 0.8 zeta 16 &
     reaction pre1 post1 react.map 1.9 taper 1.7 &
     reaction pre2 post2 react.map 1.9 taper 1.7 offset 0.843 &
     shells 1 output msevb.out 1
   fix_modify evb energy yes

.. code-block:: LAMMPS

   molecule pre pre.mol
   molecule post post.mol
   fix evb all msevb &
     coupling none &
     reaction pre post react.map 3.5 &
     shells 1 fermi_dirac 300.0 scf_topology yes &
     output msevb.out 1
   fix_modify evb energy yes

Description
"""""""""""

Perform a multi-state empirical valence bond (MSEVB) molecular dynamics
simulation.  This fix mixes the potential energy surfaces of multiple
bonding states to produce a combined surface on which the atoms propagate.
It is typically used to model proton transfer or other ion-hopping reactions
in condensed-phase systems.

Example inputs are provided in the ``examples/PACKAGES/msevb`` directory
covering the *raiteri2011*, *vuilleumier1998*, and *grimme2015* coupling
schemes, energy offsets, Fermi-Dirac smearing, and SCF topology convergence.

----------

Partitions and setup
''''''''''''''''''''

This fix is designed to be run with multiple :doc:`processor partitions
<Howto_replica>` launched via the :doc:`-partition <Run_options>`
command-line switch.  Each partition evaluates forces for one EVB bonding
state.  Using more than one partition is strongly recommended: a single
partition evaluates all states sequentially and will be significantly slower
than distributing them.  All partitions must use an identical domain
decomposition (same number of processors and box subdivision).  Positions,
box geometry, and velocities are synchronized across partitions at each step.

All reactive states are detected once per step and then distributed across
partitions in batches.  States assigned to a given partition are evaluated
in parallel with those on other partitions; if the number of states exceeds
the number of partitions, subsequent batches reuse the partitions
sequentially until all states are evaluated.  Running with more partitions
reduces the number of batches and therefore the wall time per step.

.. note::

   The pre- and post-reaction molecule templates and map file use the same
   format as :doc:`fix bond/react <fix_bond_react>`: templates are loaded
   with the :doc:`molecule <molecule>` command, and the map file specifies
   the atom-by-atom correspondence between the pre- and post-reaction
   templates and identifies the two *initiator* atoms (the transferring atom
   H and the acceptor Y).  Only basic topology changes are supported
   (bond creation and deletion, atom type and charge changes, angle, dihedral,
   and improper updates).  The more advanced bond/react features such as atom
   insertion or deletion, modification of atom IDs, and the use of atom labels
   in templates are not implemented.

----------

The EVB Hamiltonian
'''''''''''''''''''

At each timestep, this fix constructs an :math:`M \times M` effective
Hamiltonian:

.. math::

   \hat{\mathcal{H}} =
   \begin{pmatrix}
   E_{0} + \epsilon_0 & C_{0,1} & \cdots & C_{0,M-1} \\
   C_{1,0} & E_{1} + \epsilon_1 & \cdots & C_{1,M-1} \\
   \vdots  & \vdots  & \ddots & \vdots  \\
   C_{M-1,0} & C_{M-1,1} & \cdots & E_{M-1} + \epsilon_{M-1}
   \end{pmatrix}

where :math:`M` is the total number of EVB states (reference plus all
reactive states detected at this step), :math:`E_i` is the total potential
energy of state :math:`i` evaluated using that state's topology, and
:math:`\epsilon_i` is a user-specified constant energy offset (*offset*
keyword).  The off-diagonal coupling elements :math:`C_{i,j}` are
described in the Coupling section below.

Diagonalizing :math:`\hat{\mathcal{H}}` yields eigenvalues
:math:`\{\lambda_0, \ldots, \lambda_{M-1}\}` and eigenvectors
:math:`\{\mathbf{u}_0, \ldots, \mathbf{u}_{M-1}\}`.  In the default case,
the total system energy is the ground-state eigenvalue :math:`\lambda_g`
(the smallest eigenvalue), and the coefficients of the corresponding
eigenvector :math:`\mathbf{u}_g = (c_0^g, \ldots, c_{M-1}^g)` define the
amplitude of each state.  Forces are obtained from the Hellmann-Feynman
theorem:

.. math::

   \mathbf{F}_\alpha = -\frac{\partial \lambda_g}{\partial \alpha}
   = -\mathbf{u}_g^\top \frac{\partial \hat{\mathcal{H}}}{\partial \alpha}
   \mathbf{u}_g

The diagonal elements of :math:`\partial\hat{\mathcal{H}}/\partial\alpha`
are the forces from each state's force evaluation.  When the number of
states exceeds the number of partitions, some partitions evaluate more than
one state per step; when fewer states are detected than partitions, some
partitions are idle for that step.  The off-diagonal coupling forces are
computed analytically from the coupling function derivatives.

----------

Coupling
''''''''

The *coupling* keyword sets the functional form used to compute the
off-diagonal Hamiltonian elements.  A global coupling applies to all
reactions that do not specify their own; a per-reaction *coupling* keyword
inside a *reaction* block overrides it for that reaction only.  If no global
coupling is given, every *reaction* block must provide its own.

Per-reaction keywords (*coupling*, *taper*, *offset*, *shells*) may appear
in any order after the four required *reaction* arguments.

*none*

  All off-diagonal elements are zero.  The Hamiltonian is diagonal and no
  force mixing between states occurs.  Useful with *fermi_dirac* when only
  energy-based state mixing is desired, for example for polaron hopping.

*raiteri2011*

  Geometry-dependent coupling based on the asymmetric proton-transfer
  coordinate :ref:`(Raiteri2011) <msevb-Raiteri2011>`:

  .. math::

     C_{ij} = \lambda \exp\!\left(-\zeta\, Q^2\right)

  where :math:`Q = r_{HY} - r_{HX}`, with :math:`r_{HX}` the distance from
  the transferring atom H to the donor X and :math:`r_{HY}` the distance
  from H to the acceptor Y.  Forces on H, X, and Y are computed analytically.

*vuilleumier1998*

  Geometry-dependent coupling based on the donor-acceptor distance and the
  position of H relative to the midpoint of X and Y
  :ref:`(Vuilleumier1998) <msevb-Vuilleumier1998>`:

  .. math::

     C_{ij} = v_{12} \exp\!\left(-\alpha\, Q - \gamma\, q^2\right)

  where :math:`Q = r_{XY}` and :math:`q` is the distance from H to the
  midpoint of X and Y.

*grimme2015*

  Energy-dependent coupling based on the potential energy difference between
  the parent and child EVB states :ref:`(Hartke2015) <msevb-Hartke2015>`:

  .. math::

     C_{ij} = a \exp\!\left(-b\, \Delta E_{ij}^2\right)

  where :math:`\Delta E_{ij} = E_i - E_j`.

The optional *taper* keyword applies a smooth MDF damping function that
reduces the coupling to zero between the taper distance *T* and the
reaction *cutoff* distance.  The damping is based on the distance between
the initiator atoms H and Y and applies to all coupling styles, including
*grimme2015*.  This prevents discontinuities at the boundary of the
reactive zone.  The *taper* keyword may appear at the global level
(applies to all reactions without their own *taper*) or inside a
*reaction* block (per-reaction).

----------

Fermi-Dirac smearing
''''''''''''''''''''

When *fermi_dirac* is set, the system energy is represented as a
thermally weighted sum over all eigenstates rather than the ground state
alone.  Occupancies are computed self-consistently via the Fermi-Dirac
distribution:

.. math::

   f_i = \frac{1}{1 + \exp\!\left(\frac{\lambda_i - \mu}{k_B T}\right)}

where :math:`\mu` is the Fermi level solved iteratively to satisfy
:math:`\sum_i f_i = 1`.  The total energy is then
:math:`E = \sum_i f_i \lambda_i` and state amplitudes used for force
mixing are the diagonal elements of the density matrix
:math:`\rho_{ii} = \sum_j f_j |u_{ij}|^2`.

This option is useful for systems with nearly degenerate states, such
as polarons in transition-metal oxides, where the ground-state approximation
can introduce discontinuities.

----------

Multi-shell reactions
'''''''''''''''''''''

The *shells* keyword controls how many bond-hopping steps are considered
when enumerating reactive states.  At shell depth 1 (the default), only
direct donor-acceptor pairs within *cutoff* are reactive.  At depth N,
the code searches for chains of length up to N by traversing the bonded
neighbors of each detected reactive site recursively.  Each additional
shell can dramatically increase the number of EVB states.

----------

Permanent transfer
''''''''''''''''''

After each EVB evaluation, this fix checks whether the state with the
largest amplitude is no longer the reference state (state 0).  If so, it
commits the corresponding topology change permanently: bond types, atom
types, charges, angles, dihedrals, and impropers are updated, the neighbor
list is rebuilt, and the new topology is broadcast to all partitions.
A line is written to the LAMMPS log recording the reaction type and the
atom IDs involved.

When *scf_topology yes* is set, the EVB evaluation is repeated after each
permanent transfer until no further transfer is detected or *scf_max_iter*
iterations are reached.  This self-consistent loop is important when
consecutive hops are energetically favorable.

----------

Reactive atom group
'''''''''''''''''''

At initialization, this fix creates a LAMMPS group named ``<fix-ID>_atoms``
that is updated after each reactive-site detection.  The group contains
all atoms matched by the pre-reaction molecule template at any active
reactive site across all chain depths.  This group can be used in
:doc:`dump <dump>` commands to write trajectories of only the MSEVB-active
atoms:

.. code-block:: LAMMPS

   dump msevb_dump evb_atoms custom 100 msevb.lammpstrj id type xu yu zu

----------

Output file
'''''''''''

When the *output* keyword is given, a structured text file is written every
*every* steps (and always on permanent transfer).  Each record contains the
current timestep, the list of detected reactive sites, the full Hamiltonian
matrix, and the eigenvalues and amplitudes (or Fermi-Dirac occupancies if
enabled).

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  The reference topology is maintained in memory and must be
re-established from the initial data file at the start of each run.

The :doc:`fix_modify energy yes <fix_modify>` option must be set to include
the EVB energy correction (ground-state eigenvalue minus the reference
partition energy) in the system potential energy reported by thermo.

This fix stores a global scalar (the EVB energy correction,
:math:`\lambda_g - E_0`) and a global vector of length equal to the number
of partitions containing the potential energy of each partition.  Both the
scalar and vector are in :doc:`energy units <units>` and "extensive".
See the :doc:`Howto output <Howto_output>` page for how to access these
values.

The virial contribution from the coupling forces is tallied and included
when :doc:`fix_modify virial yes <fix_modify>` is set.

This fix is invoked during :doc:`energy minimization <minimize>`.  The EVB
force mixing and permanent transfer logic runs at each minimization step.
Note that permanent topology changes during minimization trigger a neighbor
list rebuild and effectively restart the local minimization from the new
topology.

Restrictions
""""""""""""

This fix is part of the MSEVB package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

Exactly one instance of this fix may be active at a time.

All partitions must have an identical domain decomposition (same number
of processors and the same spatial subdivision).

This fix is not compatible with :doc:`fix balance <fix_balance>` or
the *respa* run style.

Related commands
""""""""""""""""

:doc:`fix bond/react <fix_bond_react>` command,
:doc:`fix alchemy <fix_alchemy>` command,
:doc:`molecule <molecule>` command,
:doc:`-partition command-line switch <Run_options>`

Default
"""""""

The default values are *shells* = 1, *scf_topology* = no,
*scf_max_iter* = 10, and *fermi_dirac* disabled.
No global coupling is set; every *reaction* block must provide its own
if no global *coupling* keyword is given.

----------

.. _msevb-Raiteri2011:

**(Raiteri2011)** Raiteri P, Gale JD, Bussi G. Reactive force field simulation of proton diffusion in BaZrO3 using an empirical valence bond approach. Journal of Physics: Condensed Matter. 2011 Aug 24;23(33):334213.

.. _msevb-Vuilleumier1998:

**(Vuilleumier1998)** Vuilleumier R, Borgis D. An extended empirical valence bond model for describing proton transfer in H+ (H2O) n clusters and liquid water. Chemical physics letters. 1998 Feb 20;284(1-2):71-7.

.. _msevb-Hartke2015:

**(Hartke2015)** Hartke B, Grimme S. Reactive force fields made simple. Physical Chemistry Chemical Physics. 2015;17(26):16715-8.
