.. index:: fix_modify AtC reset_atomic_reference_positions

fix_modify AtC reset_atomic_reference_positions command
=======================================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> reset_atomic_reference_positions

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* reset_atomic_reference_positions = name of the AtC sub-command


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC reset_atomic_reference_positions

Description
"""""""""""

Resets the atomic positions ATC uses to perform point to field
operations. In can be used to use perfect lattice sites in ATC but a
thermalized or deformed lattice in LAMMPS.


Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

None

