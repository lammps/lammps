.. index:: fix_modify AtC add_molecule

fix_modify AtC add_molecule command
===================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> add_molecule <small|large> <tag> <group-ID>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* add_molecule = name of the AtC sub-command
* *small* or *large* = can be *small* if molecule size < cutoff radius, must be *large* otherwise
* tag = tag for tracking a molecule
* *group-ID* = LAMMPS defined group-ID

Examples
""""""""

.. code-block:: LAMMPS

   group WATERGROUP type 1 2
   fix_modify AtC add_molecule small water WATERGROUP

Description
"""""""""""

Associates a tag with all molecules corresponding to a specified group.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC add_species <atc_add_species>`
- :doc:`fix_modify AtC remove_species <atc_remove_species>`
- :doc:`fix_modify AtC remove_molecule <atc_remove_molecule>`

Default
"""""""

None.
