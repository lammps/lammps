.. index:: fix_modify AtC unfix

fix_modify AtC initial command
==============================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> unfix <field> <nodeset>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* unfix = name of the AtC sub-command
* field = field kind name valid for type of physics: *temperature* or *electron_temperature*
* nodeset = name of set of nodes to apply boundary condition

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC unfix temperature groupNAME

Description
"""""""""""

Removes constraint on field values for specified nodes.

Restrictions
""""""""""""

The keyword *all* is reserved and thus not available as nodeset name.

Related AtC commands
""""""""""""""""""""

:doc:`fix_modify AtC fix <atc_fix>`

Default
"""""""

None.
