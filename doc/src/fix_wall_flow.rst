.. index:: fix wall/flow
.. index:: fix wall/flow/kk

fix wall/flow command
=====================

Accelerator Variants: *wall/flow/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID wall/flow axis vflow T seed N coords ... keyword value

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/flow = style name of this fix command
* axis = flow axis (*x*, *y*, or *z*)
* vflow = generated flow velocity in *axis* direction (velocity units)
* T = flow temperature (temperature units)
* seed = random seed for stochasticity (positive integer)
* N = number of walls
* coords = list of N wall positions along the *axis* direction in ascending order (distance units)
* zero or more keyword/value pairs may be appended
* keyword = *units*

  .. parsed-literal::

       *units* value = *lattice* or *box*
         *lattice* = wall positions are defined in lattice units
         *box* = the wall positions are defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all wall/flow x 0.4 1.5 593894 4 2.0 4.0 6.0 8.0

Description
"""""""""""

.. versionadded:: TBD

This fix implements flow boundary conditions (FBC) introduced in
:ref:`(Pavlov1) <fbc-Pavlov1>` and :ref:`(Pavlov2) <fbc-Pavlov2>`.
The goal is to generate a stationary flow with a shifted Maxwell
velocity distribution:

.. math::

   f_a(v_a) \propto \exp{\left(-\frac{m (v_a-v_{\text{flow}})^2}{2 kB T}\right)}

where :math:`v_a` is the component of velocity along the specified
*axis* argument (a = x,y,z), :math:`v_{\text{flow}}` is the flow
velocity specified as the *vflow* argument, *T* is the specified flow
temperature, *m* is the particle mass, and *kB* is the Boltzmann
constant.

This is achieved by defining a series of *N* transparent walls along
the flow *axis* direction.  Each wall is at the specified position
listed in the *coords* argument.  Note that an additional transparent
wall is defined by the code at the boundary of the (periodic)
simulation domain in the *axis* direction.  So there are effectively
N+1 walls.

Each time a particle in the specified group passes through one of the
transparent walls, its velocity is re-assigned.  Particles not in the
group do not interact with the wall. This can be used, for example, to
add obstacles composed of atoms, or to simulate a solution of complex
molecules in a one-atom liquid (note that the fix has been tested for
one-atom systems only).

Conceptually, the velocity re-assignment represents creation of a new
particle within the system with simultaneous removal of the particle
which passed through the wall.  The velocity components in directions
parallel to the wall are re-assigned according to the standard Maxwell
velocity distribution for the specified temperature *T*.  The velocity
component perpendicular to the wall is re-assigned according to the
shifted Maxwell distribution defined above:

.. math::

   f_{\text{a generated}}(v_a) \propto v_a f_a(v_a)

It can be shown that for an ideal-gas scenario this procedure makes
the velocity distribution of particles between walls exactly as
desired.

Since in most cases simulated systems are not an ideal gas, multiple
walls can be defined, since a single wall may not be sufficient for
maintaining a stationary flow without "congestion" which can manifest
itself as regions in the flow with increased particle density located
upstream from static obstacles.

For the same reason, the actual temperature and velocity of the
generated flow may differ from what is requested.  The degree of
discrepancy is determined by how different from an ideal gas the
simulated system is.  Therefore, a calibration procedure may be
required for such a system as described in :ref:`(Pavlov)
<fbc-Pavlov2>`.

Note that the interactions between particles on different sides of a
transparent wall are not disabled or neglected.  Likewise particle
positions are not altered by the velocity reassignment.  This removes
the need to modify the force field to work correctly in cases when a
particle is close to a wall.

For example, if particle positions were uniformly redistributed across
the surface of a wall, two particles could end up too close to each
other, potentially causing the simulation to explode.  However due to
this compromise, some collective phenomena such as regions with
increased/decreased density or collective movements are not fully
removed when particles cross a wall.  This unwanted consequence can
also be potentially mitigated by using more multiple walls.

.. note::

  When the specified flow has a high velocity, a lost atoms error can
  occur (see :doc:`error messages <Errors_messages>`).  If this
  happens, you should ensure the checks for neighbor list rebuilds,
  set via the :doc:`neigh_modify <neigh_modify>` command, are as
  conservative as possible (every timestep if needed).  Those are the
  default settings.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to
this fix.

No global or per-atom quantities are stored by this fix for access by
various :doc:`output commands <Howto_output>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

Fix *wall_flow* is part of the EXTRA-FIX package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Flow boundary conditions should not be used with rigid bodies such as
those defined by a "fix rigid" command.

This fix can only be used with periodic boundary conditions along the
flow axis. The size of the box in this direction must not change. Also,
the fix is designed to work only in an orthogonal simulation box.

Related commands
""""""""""""""""

:doc:`fix wall/reflect <fix_wall>` command

Default
"""""""

The default for the units keyword is lattice.

----------

.. _fbc-Pavlov1:

**(Pavlov1)** Pavlov, Kolotinskii, Stegailov, "GPU-Based Molecular Dynamics of Turbulent Liquid Flows with OpenMM", Proceedings of PPAM-2022, LNCS (Springer), vol. 13826, pp. 346-358 (2023)

.. _fbc-Pavlov2:

**(Pavlov2)** Pavlov, Galigerov, Kolotinskii, Nikolskiy, Stegailov, "GPU-based Molecular Dynamics of Fluid Flows: Reaching for Turbulence", Int. J. High Perf. Comp. Appl., (2024)
