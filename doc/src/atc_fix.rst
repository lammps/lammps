.. index:: fix_modify AtC fix

fix_modify AtC fix command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> fix <field> <nodeset> <constant|function>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* fix = name of the AtC sub-command
* field = field kind name valid for type of physics: *temperature* or *electron_temperature*
* nodeset = name of set of nodes to apply boundary condition
* *constant* or *function* = value or name of function followed by its parameters

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC fix temperature groupNAME 10.
   fix_modify AtC fix temperature groupNAME 0 0 0 10.0 0 0 1.0

Description
"""""""""""

Creates a constraint on the values of the specified field at specified nodes.

Restrictions
""""""""""""

The keyword *all* is reserved and thus not available as nodeset name.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC unfix <atc_unfix>`

Default
"""""""

None.
