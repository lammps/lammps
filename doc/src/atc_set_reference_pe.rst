.. index:: fix_modify AtC set reference_potential_energy

fix_modify AtC set reference_potential_energy command
=====================================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> set reference_potential_energy [<value|filename>]

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* set reference_potential_energy = name of the AtC sub-command
* value = optional user specified zero point for PE in native LAMMPS energy units
* filename = optional user specified string for file of nodal PE values to be read-in

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC set reference_potential_energy
   fix_modify AtC set reference_potential_energy -0.05
   fix_modify AtC set reference_potential_energy myPEvalues

Description
"""""""""""

Used to set various quantities for the post-processing algorithms. It
sets the zero point for the potential energy density using the value
provided for all nodes, or from the current configuration of the lattice
if no value is provided, or values provided within the specified
filename.

Restrictions
""""""""""""

Must be used with :doc:`fix atc hardy <fix_atc>` or :doc:`fix atc field <fix_atc>`.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

Defaults to the LAMMPS zero point i.e. isolated atoms.

