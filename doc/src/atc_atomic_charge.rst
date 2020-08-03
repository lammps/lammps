.. index:: fix_modify AtC atomic_charge

fix_modify AtC atomic_charge command
====================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> <include|omit> atomic_charge

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* *include* or *omit* = switch to activate/deactivate inclusion of intrinsic atomic charge in ATC
* atomic_charge = name of the AtC sub-command

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC include atomic_charge

Description
"""""""""""

Determines whether AtC tracks the total charge as a finite element
field.

Restrictions
""""""""""""

Required for: *electrostatics*

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

If the atom charge is defined, default is on, otherwise default is off.
