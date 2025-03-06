.. index:: fix_modify AtC add_species

fix_modify AtC add_species command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> add_species <tag> <group|type> <ID>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* add_species = name of the AtC sub-command
* tag = tag for tracking a species
* *group* or *type* = LAMMPS defined group or type of atoms
* ID = name of group or type number

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC add_species gold type 1
   group GOLDGROUP type 1
   fix_modify AtC add_species gold group GOLDGROUP

Description
"""""""""""

Associates a tag with all atoms of a specified type or within a specified group.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC add_molecule <atc_add_molecule>`
- :doc:`fix_modify AtC remove_species <atc_remove_species>`
- :doc:`fix_modify AtC remove_molecule <atc_remove_molecule>`

Default
"""""""

None.
