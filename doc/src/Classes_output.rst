Output Class
************

The Output class manages all forms of output from LAMMPS simulations, including
thermodynamic data, dump files, and restart files.  It coordinates the timing
and execution of output operations during simulation runs.

Key responsibilities include:

- Managing thermodynamic output (thermo) and its frequency
- Managing multiple dump files with different output frequencies
- Managing restart file output with single or double file modes
- Coordinating output timing to determine the next output timestep
- Supporting variable-based output frequencies
- Providing a factory for creating dump style instances

Output timing can be based on timestep number or simulation time.  The class
tracks the next timestep for each type of output and provides methods to force
immediate output of dumps or restart files.

Thermodynamic output is controlled by the :doc:`thermo <thermo>`,
:doc:`thermo_style <thermo_style>`, and :doc:`thermo_modify <thermo_modify>`
commands.  Dump files are created with the :doc:`dump <dump>` command and
restart files with the :doc:`restart <restart>` command.

.. doxygenclass:: LAMMPS_NS::Output
   :project: progguide
   :members:
