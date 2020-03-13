.. index:: fix_modify AtC time_integration

fix_modify AtC time_integration command
=============================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> time_integration <descriptor>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* time\_integration = name of the AtC sub-command
* descriptor =  *gear* or *fractional_step* or *verlet*


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify atc time_integration fractional_step


Description
"""""""""""

Command to select the thermal or momentum time integration.

---------

Options for thermal time integration:

*gear*
  atomic velocity update with 2nd order Verlet, nodal temperature update
  with 3rd or 4th order Gear, thermostats based on controlling power

*fractional\_step*
  atomic velocity update with 2nd order Verlet, mixed nodal temperature
  update, 3/4 Gear for continuum and 2 Verlet for atomic contributions,
  thermostats based on controlling discrete energy changes

---------

Options for momentum time integration:

*verlet*
  atomic velocity update with 2nd order Verlet, nodal temperature update
  with 2nd order Verlet, kinetostats based on controlling force

*fractional\_step*
  atomic velocity update with 2nd order Verlet, mixed nodal momentum
  update, 2nd order Verlet for continuum and exact 2nd order Verlet for
  atomic contributions, kinetostats based on controlling discrete
  momentum changes

*gear*
  atomic velocity update with 2nd order Verlet, nodal temperature update
  with 3rd or 4th order Gear, kinetostats based on controlling power.

---------

Restrictions
""""""""""""

None.

Related commands
""""""""""""""""

:doc:`fix atc <fix_atc>`

Default
"""""""

None.
