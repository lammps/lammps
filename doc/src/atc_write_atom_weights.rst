.. index:: fix_modify AtC write_atom_weights

fix_modify AtC write_atom_weights command
=========================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> write_atom_weights <filename> <frequency>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* write_atom_weights = name of the AtC sub-command
* filename = name of file that atomic weights are written to
* frequency = how often writes will occur

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC write_atom_weights atm_wt_file.txt 10

Description
"""""""""""

Command for writing the values of atomic weights to a specified file.


Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

None

