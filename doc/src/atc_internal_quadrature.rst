.. index:: fix_modify AtC internal_quadrature

fix_modify AtC internal_quadrature command
==========================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> internal_quadrature <on|off> [region]

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* internal_quadrature = name of the AtC sub-command
* on or off = select whether internal quadrature is enabled or not
* region = treat finite elements as within MD region (optional)


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC internal_quadrature off

Description
"""""""""""

Command to use or not use atomic quadrature on internal elements fully
filled with atoms. By turning the internal quadrature off these elements
do not contribute to the governing PDE and the fields at the internal
nodes follow the weighted averages of the atomic data.

Optional region tag specifies which finite element nodes will be treated
as being within the MD region. This option is only valid with
internal_quadrature off.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

on.
