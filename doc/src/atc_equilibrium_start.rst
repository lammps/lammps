.. index:: fix_modify AtC equilibrium_start

fix_modify AtC equilibrium_start command
========================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> equilibrium_start <on|off>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* equilibrium_start = name of the AtC sub-command
* *exponential* or *step* or *no_filter* = select type of filter

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC equilibrium_start on

Description
"""""""""""

Starts filtered calculations assuming they start in equilibrium,
i.e. perfect finite element force balance.

Restrictions
""""""""""""

Only for use with these specific transfers: thermal, two_temperature

Related AtC commands
""""""""""""""""""""
- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC filter <atc_time_filter>`
- :doc:`fix_modify AtC filter scale <atc_filter_scale>`

Default
"""""""

None.
