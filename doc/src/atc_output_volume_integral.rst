.. index:: fix_modify AtC output volume_integral

fix_modify AtC output volume_integral command
=============================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> output volume_integral <elementset_name> <field>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* output volume_integral = name of the AtC sub-command
* elementset_name= name of elementset to be integrated over
* fieldname = name of field to integrate


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC output volume_integral eset1 mass_density


Description
"""""""""""

Performs volume integration of specified field over elementset and
outputs resulting variable values to GLOBALS file.

Restrictions
""""""""""""

None.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix atc command <fix_atc>`

Default
"""""""

None.
