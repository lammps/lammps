.. index:: fix wall/flow
.. index:: fix wall/flow/kk

fix wall/flow command
=====================

Accelerator Variants: *wall/flow/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID wall/flow ax vf T seed N coords ... keyword value

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/flow = style name of this fix command
* ax = flow axis (*x*, *y*, or *z* character)
* vf = *ax* component of generated flow velocity
* T = flow temperature (temperature units)
* seed = random seed for stochasticity (positive integer)
* N = number of walls (positive integer)
* coords = set of N wall coordinates (box units) along *ax* axis arranged in ascending order. Note that an additional implicit wall is introduced at the boundary of the simulation domain, so the resulting system always has N+1 walls.

* zero or more keyword/value pairs may be appended
* keyword = *units*

  .. parsed-literal::

       *units* value = *lattice* or *box*
         *lattice* = the wall positions are defined in lattice units
         *box* = the wall positions are defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 g_flow wall/flow x ${VFLOW} ${TEMP} 123 ${nwall} ${w1} ${w2} ${w3} ${w4}
   fix 2 all wall/flow 0.4 0.2 3 1 400

Description
"""""""""""

This fix implements flow boundary conditions (FBC) introduced in :ref:`(Pavlov1) <fbc-Pavlov1>` and :ref:`(Pavlov2) <fbc-Pavlov2>`.
The goal is to generate a stationary flow with a shifted Maxwell velocity distribution:

.. math::

   f_z(v_z) \propto \exp{\left(-\frac{m (v_z-v_{\text{flow}})^2}{2 k T}\right)}

This is achieved by reassigning the velocity of each particle that passes a wall.
Such reassigning represents an emission of a new particle into the system with
simultaneous removal of a particle with the same position.
The velocity components parallel to the wall are re-assigned according
to the Maxwell velocity distribution. The perpendicular component is assigned
according to the following velocity distribution:

.. math::

   f_{\text{z generated}}(v_z) \propto v_z f_z(v_z)

It can be shown that in an ideal-gas scenario this makes the velocity
distribution of particles between walls exactly as desired.

Since in most cases simulated systems are not ideal gas,
the need for multiple walls might arise, as a single wall may not be
sufficient for maintaining a stationary flow without congestion
manifesting as areas with increased density located upstream from static obstacles.

For the same reason, the actual temperature and velocity of the generated
flow may differ from ones requested. The degree of such discrepancy is determined
by how different from the ideal gas the simulated system is. Therefore, a calibration procedure is required for each system as described in :ref:`(Pavlov) <fbc-Pavlov2>`.

The interactions between particles on different sides of a wall are not disabled or neglected and the
particle positions are not affected by the velocity reassignment.
This removes the need to modify the force field to work correctly in cases when a particle is close
to a wall (for example, if particle positions were uniformly redistributed across the surface of the wall,
two particles could end up too close to each other, potentially causing the simulation to explode).
However due to this compromise, some collective phenomena such as areas with increased/decreased density
or collective movements are not fully removed when particles cross a wall.
This unwanted consequence can also be potentially mitigated by using more than one wall.


----------

Note that when high flow velocity is reached, a lost atoms error may
occur (see :doc:`error messages <Errors_messages>`).
If this message appears when using this fix, you can, for example, reduce the frequency of the
neighbor list rebuild via :doc:`neigh_modify <neigh_modify>` command.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

None of the :doc:`fix_modify <fix_modify>` options are relevant to
this fix.

No global or per-atom quantities are stored by this fix for access by
various :doc:`output commands <Howto_output>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

Flow boundary conditions should not be used with rigid bodies such as those
defined by a "fix rigid" command.

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
