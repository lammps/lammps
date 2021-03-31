.. index:: fix_modify AtC time_integration

fix_modify AtC time_integration command
=======================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> time_integration <descriptor>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* time_integration = name of the AtC sub-command
* descriptor =  *gear* or *fractional_step* or *verlet*


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC time_integration fractional_step


Description
"""""""""""

Command to select the thermal or momentum time integration.

---------

Options for thermal time integration:

*gear*
  atomic velocity update with second order Verlet, nodal temperature update
  with third or fourth order Gear, thermostats based on controlling power

*fractional_step*
  atomic velocity update with second order Verlet, mixed nodal temperature
  update, 3/4 Gear for continuum and 2 Verlet for atomic contributions,
  thermostats based on controlling discrete energy changes

---------

Options for momentum time integration:

*verlet*
  atomic velocity update with second order Verlet, nodal temperature update
  with second order Verlet, kinetostats based on controlling force

*fractional_step*
  atomic velocity update with second order Verlet, mixed nodal momentum
  update, second order Verlet for continuum and exact second order Verlet for
  atomic contributions, kinetostats based on controlling discrete
  momentum changes

*gear*
  atomic velocity update with second order Verlet, nodal temperature update
  with third or fourth order Gear, kinetostats based on controlling power.

---------

Restrictions
""""""""""""

None.


Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

None.
