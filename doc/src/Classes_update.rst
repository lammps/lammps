Update Class
************

The Update class manages time integration (molecular dynamics) and
energy minimization in LAMMPS.  It holds the simulation timestep size
(``dt``, which can be (re-)set with the :doc:`timestep command
<timestep>`) and tracks the current timestep number (``ntimestep``,
which can be reset with the :doc:`reset_timestep command
<reset_timestep>`).  That simulation set is recorded using the
``whichflag`` member (1 for MD, 2 for minimization, and 0 for neither).

Key responsibilities include:

- Managing the timestep size and current step number
- Creating and managing time integrators (e.g., Verlet, rRESPA)
- Creating and managing energy minimizers (e.g., CG, FIRE, SD)
- Tracking simulation time and run boundaries
- Managing energy and virial tallying flags
- Handling unit style settings that affect the timestep

The class maintains factory maps for dynamically creating integrator and
minimizer styles based on user input from the :doc:`run_style <run_style>`
and :doc:`min_style <min_style>` commands.

.. doxygenclass:: LAMMPS_NS::Update
   :project: progguide
   :members:
