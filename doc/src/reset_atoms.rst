.. index:: reset_atoms

reset_atoms command
===================

Syntax
""""""

.. code-block:: LAMMPS

   reset_atoms sub-command keyword values ...

   * sub-command = *id* or *image* or *mol*
   * zero or more keyword/value pairs may be appended depending on sub-command

Examples
""""""""

.. code-block:: LAMMPS

   reset_atoms id
   reset_atoms image all
   reset_atoms mol all

Description
"""""""""""

.. versionadded:: TBD

The *reset_atoms* command provides a number of sub-commands that reset
selected atom properties like atom IDs, molecule IDs, or image flags
according to selected algorithms.  Those are often useful when the
simulated system has undergone some significant modifications like
adding or removing atoms or molecules, joining data files, changing
bonds, or diffusion.  Follow the links listed below to see the
documentation for individual sub-commands.

- :doc:`reset_atoms_id`
- :doc:`reset_atoms_image`
- :doc:`reset_atoms_mol`


Defaults
""""""""

none
