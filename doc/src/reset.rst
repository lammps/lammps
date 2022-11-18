.. index:: reset

reset command
=============

Syntax
""""""

.. code-block:: LAMMPS

   reset sub-command keyword values ...

   * sub-command = *atom_ids* or *image_flags* or *mol_ids*
   * zero or more keyword/value pairs may be appended depending on sub-command

Examples
""""""""

.. code-block:: LAMMPS

   reset atom_ids
   reset image_flags all
   reset mol_ids all

Description
"""""""""""

.. versionadded:: TBD

The *reset* command provides a number of sub-commands that reset
selected atom properties like atom IDs, molecule IDs, or image flags
according to selected algorithms.  Those are often useful when the
simulated system has undergone some significant modifications like
adding or removing atoms or molecules, joining data files, changing
bonds, or diffusion.  Follow the links listed below to see the
documentation for individual sub-commands.

- :doc:`reset_atom_ids`
- :doc:`reset_image_flags`
- :doc:`reset_mol_ids`


Defaults
""""""""

none
