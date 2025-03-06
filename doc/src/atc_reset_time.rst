.. index:: fix_modify AtC reset_time

fix_modify AtC reset_time command
=================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> reset_time <value>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* reset_time = name of the AtC sub-command
* value = new time value

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC reset_time 0.0

Description
"""""""""""

Resets the simulation time counter.


Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

None

