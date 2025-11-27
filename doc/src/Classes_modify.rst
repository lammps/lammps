LAMMPS Modify Class
*******************

The Modify class manages all :doc:`Fix <fix>` and :doc:`Compute <compute>`
instances in a LAMMPS simulation.  Fixes modify the simulation state
(e.g., applying forces, constraints, or thermostats), while Computes
calculate properties from the simulation (e.g., temperature, pressure,
or per-atom quantities).

Key responsibilities include:

- Creating, storing, and deleting Fix and Compute instances
- Calling fixes at appropriate points during the timestep (initial_integrate,
  post_force, final_integrate, etc.)
- Managing fix execution order via bitmask flags
- Handling restart file read/write for fixes
- Managing energy contributions from fixes for thermodynamics
- Coordinating minimization callbacks for fixes

The class maintains lists of fixes to be called at each stage of the
timestep, organized by callback type.  It also provides factory maps for
dynamically creating new fix and compute styles.

.. doxygenclass:: LAMMPS_NS::Modify
   :project: progguide
   :members:
