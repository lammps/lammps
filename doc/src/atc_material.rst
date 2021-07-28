.. index:: fix_modify AtC material

fix_modify AtC material command
===============================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> material <elementset_name> <material_id>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* material = name of the AtC sub-command
* elementset_name = name of the elementset
* material_id = ID of the material

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC material gap_region 1

Description
"""""""""""

Sets the material model in *elementset_name* to be of type *material_id*\ .

Restrictions
""""""""""""

The element set must already be created and the material must be
specified in the material file given the the atc fix on construction

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

All elements default to the first material in the material file.
