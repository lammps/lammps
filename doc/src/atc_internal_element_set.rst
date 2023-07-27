.. index:: fix_modify AtC internal_element_set

fix_modify AtC internal_element_set command
===========================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> internal_element_set <element_set_name>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* internal_element_set = name of the AtC sub-command
* element_set_name = name of element set defining internal region, or *off*

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC internal_element_set myElementSet
   fix_modify AtC internal_element_set off

Description
"""""""""""

Enables AtC to base the region for internal atoms to be an element
set. If no ghost atoms are used, all the AtC atoms must be constrained
to remain in this element set by the user, e.g., with walls. If boundary
atoms are used in conjunction with Eulerian atom maps AtC will partition
all atoms of a boundary or internal type to be of type internal if they
are in the internal region or to be of type boundary otherwise.

Restrictions
""""""""""""

If boundary atoms are used in conjunction with Eulerian atom maps, the
Eulerian reset frequency must be an integer multiple of the Lammps
reneighbor frequency.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC atom_element_map <atc_atom_element_map>`
- :doc:`fix_modify AtC boundary type <atc_boundary_type>`

Default
"""""""

*off*
