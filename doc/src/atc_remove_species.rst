.. index:: fix_modify AtC remove_species

fix_modify AtC remove_species command
=====================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> remove_species <tag>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* remove_species = name of the AtC sub-command
* tag = tag for tracking a species

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC remove_species gold

Description
"""""""""""

Removes tag designated for tracking a specified species.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC add_species <atc_add_species>`
- :doc:`fix_modify AtC add_molecule <atc_add_molecule>`
- :doc:`fix_modify AtC remove_molecule <atc_remove_molecule>`

Default
"""""""

None.
