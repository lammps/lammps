.. index:: set_time

set_time command
================

Syntax
""""""

.. parsed-literal::

   set_time t

* t = current simulation time (time units)

Examples
""""""""

.. code-block:: LAMMPS

   time 0.0
   time 10.5

Description
"""""""""""

Set the current accumulated simulation time for subsequent molecular
dynamics simulations.  See the :doc:`units <units>` command for the time
units associated with each choice of units that LAMMPS supports.


Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`reset_timestep <reset_timestep>`, :doc:`timestep <timestep>`,
:doc:`fix dt/reset <fix_dt_reset>`, :doc:`units <units>`

Default
"""""""

0.0 at the beginning of the first run.
