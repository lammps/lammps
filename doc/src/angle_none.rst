.. index:: angle_style none

angle_style none command
========================

Syntax
""""""

.. code-block:: LAMMPS

   angle_style none

Examples
""""""""

.. code-block:: LAMMPS

   angle_style none

Description
"""""""""""

Using an angle style of none means angle forces and energies are not
computed, even if triplets of angle atoms were listed in the data file
read by the :doc:`read_data <read_data>` command.

See the :doc:`angle_style zero <angle_zero>` command for a way to
calculate angle statistics, but compute no angle interactions.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`angle_style zero <angle_zero>`

**Default:** none
