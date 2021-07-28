.. index:: timestep

timestep command
================

Syntax
""""""

.. parsed-literal::

   timestep dt

* dt = timestep size (time units)

Examples
""""""""

.. code-block:: LAMMPS

   timestep 2.0
   timestep 0.003

Description
"""""""""""

Set the timestep size for subsequent molecular dynamics simulations.
See the :doc:`units <units>` command for the time units associated with
each choice of units that LAMMPS supports.

The default value for the timestep size also depends on the choice of
units for the simulation; see the default values below.

When the :doc:`run style <run_style>` is *respa*, dt is the timestep for
the outer loop (largest) timestep.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix dt/reset <fix_dt_reset>`, :doc:`run <run>`,
:doc:`run_style <run_style>` respa, :doc:`units <units>`

Default
"""""""

+--------------------------------+---------------+-----------------------+
| choice of :doc:`units <units>` | time units    | default timestep size |
+--------------------------------+---------------+-----------------------+
| lj                             | :math:`\tau`  | 0.005 :math:`\tau`    |
+--------------------------------+---------------+-----------------------+
| real                           | fs            | 1.0 fs                |
+--------------------------------+---------------+-----------------------+
| metal                          | ps            | 0.001 ps              |
+--------------------------------+---------------+-----------------------+
| si                             | s             | 1.0e-8 s (10 ns)      |
+--------------------------------+---------------+-----------------------+
| cgs                            | s             | 1.0e-8 s (10 ns)      |
+--------------------------------+---------------+-----------------------+
| electron                       | fs            | 0.001 fs              |
+--------------------------------+---------------+-----------------------+
| micro                          | :math:`\mu`\ s| 2.0 :math:`\mu`\ s    |
+--------------------------------+---------------+-----------------------+
| nano                           | ns            | 0.00045 ns            |
+--------------------------------+---------------+-----------------------+
