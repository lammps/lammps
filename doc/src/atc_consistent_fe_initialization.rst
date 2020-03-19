.. index:: fix_modify AtC consistent_fe_initialization

fix_modify AtC consistent_fe_initialization command
===================================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> consistent_fe_initialization <on|off>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* consistent_fe_initialization = name of the AtC sub-command
* *on* or *off* = switch to activiate/deactiviate the initial setting of FE intrinsic field to match the projected MD field

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC consistent_fe_initialization on

Description
"""""""""""

Determines whether AtC initializes FE intrinsic fields (e.g.,
temperature) to match the projected MD values. This is particularly
useful for fully overlapping simulations.

Restrictions
""""""""""""

Can be used with: *thermal*, *two_temperature*.
Cannot be used with time filtering on.
Does not include boundary nodes.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

Default is *off*
