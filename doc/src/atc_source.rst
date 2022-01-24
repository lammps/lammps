.. index:: fix_modify AtC source

fix_modify AtC source command
==============================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> source <field> <element_set> <value|function>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* source = name of the AtC sub-command
* field = field kind name valid for type of physics: *temperature* or *electron_temperature*
* element_set = name of set of elements
* *value* or *function* = value or name of function followed by its parameters

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC source temperature middle temporal_ramp 10.0 0.0

Description
"""""""""""

Add domain sources to the mesh. The units are consistent with LAMMPS's
units for mass, length and time and are defined by the PDE being solved,
e.g. for thermal transfer the balance equation is for energy and source
is energy per time.

Restrictions
""""""""""""

The keyword *all* is reserved and thus not available as element_set name.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC remove_source <atc_remove_source>`

Default
"""""""

None.
