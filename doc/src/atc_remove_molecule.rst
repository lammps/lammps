.. index:: fix_modify AtC remove_molecule

fix_modify AtC remove_molecule command
======================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> remove_molecule <tag>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* remove_molecule = name of the AtC sub-command
* tag = tag for tracking a molecule

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC remove_molecule water

Description
"""""""""""

Removes tag designated for tracking a specified set of molecules.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC add_species <atc_add_species>`
- :doc:`fix_modify AtC add_molecule <atc_add_molecule>`
- :doc:`fix_modify AtC remove_species <atc_remove_species>`

Default
"""""""

None.
