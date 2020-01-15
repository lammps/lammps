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


.. parsed-literal::

   timestep 2.0
   timestep 0.003

Description
"""""""""""

Set the timestep size for subsequent molecular dynamics simulations.
See the :doc:`units <units>` command for the time units associated with
each choice of units that LAMMPS supports.

The default value for the timestep size also depends on the choice of
units for the simulation; see the default values below.

When the :doc:`run style <run_style>` is *respa*\ , dt is the timestep for
the outer loop (largest) timestep.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix dt/reset <fix_dt_reset>`, :doc:`run <run>`,
:doc:`run\_style <run_style>` respa, :doc:`units <units>`

Default
"""""""

+--------------------------------+------------+-----------------------+
| choice of :doc:`units <units>` | time units | default timestep size |
+--------------------------------+------------+-----------------------+
| lj                             | tau        | 0.005 tau             |
+--------------------------------+------------+-----------------------+
| real                           | fmsec      | 1.0 fmsec             |
+--------------------------------+------------+-----------------------+
| metal                          | psec       | 0.001 psec            |
+--------------------------------+------------+-----------------------+
| si                             | sec        | 1.0e-8 sec (10 nsec)  |
+--------------------------------+------------+-----------------------+
| cgs                            | sec        | 1.0e-8 sec (10 nsec)  |
+--------------------------------+------------+-----------------------+
| electron                       | fmsec      | 0.001 fmsec           |
+--------------------------------+------------+-----------------------+
| micro                          | usec       | 2.0 usec              |
+--------------------------------+------------+-----------------------+
| nano                           | nsec       | 0.00045 nsec          |
+--------------------------------+------------+-----------------------+


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
