.. index:: fix_modify AtC boundary type

fix_modify AtC boundary type command
=============================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> boundary type <atom-type-id>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* boundary type = name of the AtC sub-command
* atom-type-id =  type id for atoms that represent a fictitious boundary internal to the FE mesh


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC boundary type ghost_atoms

Description
"""""""""""

Command to define the atoms that represent the ficticious boundary internal to the FE mesh. For fully overlapped MD/FE domains with periodic boundary conditions no boundary atoms should be defined.

Restrictions
""""""""""""

None.

Default
"""""""

None.
