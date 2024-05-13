.. index:: fix_modify AtC output boundary_integral

fix_modify AtC output boundary_integral command
===============================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> output boundary_integral <fieldname> faceset [name]

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* output boundary_integral = name of the AtC sub-command
* fieldname = name of hardy field
* faceset = required keyword
* name= name of faceset


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC output boundary_integral stress faceset loop1


Description
"""""""""""

Calculates a surface integral of the given field dotted with the outward
normal of the faces and puts output in the "GLOBALS" file.

Restrictions
""""""""""""

Must be used with the hardy/field type of :doc:`fix atc <fix_atc>`

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

None.
