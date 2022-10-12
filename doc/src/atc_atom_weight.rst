.. index:: fix_modify AtC atom_weight

fix_modify AtC atom_weight command
==================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> atom_weight <method> <args>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* atom_weight = name of the AtC sub-command
* method = *constant* or *lattice* or *element* or *region* or *group* or *read_in*

  - *constant* <group-ID> <value>: atoms in specified group are assigned the constant value given
  - *lattice*\ : volume per atom for specified lattice type (e.g. fcc) and parameter
  - *element*\ : element volume divided among atoms within element
  - *region*\ : volume per atom determined based on the atom count in the MD regions and their volumes. Note: meaningful only if atoms completely fill all the regions.
  - *group*\ : volume per atom determined based on the atom count in a group and its volume
  - *node*\ : (undocumented)
  - *node_element*\ : (undocumented)
  - *read_in*\ <filename>: list of values for atoms are read-in from specified file

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC atom_weight constant myatoms 11.8
   fix_modify AtC atom_weight lattice
   fix_modify AtC atom_weight read-in atm_wt_file.txt

Description
"""""""""""

Command for assigning the value of atomic weights used for atomic
integration in atom-continuum coupled simulations.


Restrictions
""""""""""""

The use of the lattice option requires a lattice type and parameter is already specified.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

*lattice*

