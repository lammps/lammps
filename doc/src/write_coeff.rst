.. index:: write_coeff

write_coeff command
===================

Syntax
""""""

.. code-block:: LAMMPS

   write_coeff file

* file = name of data file to write out

Examples
""""""""

.. code-block:: LAMMPS

   write_coeff polymer.coeff

Description
"""""""""""

Write a text format file with the currently defined force field
coefficients in a way, that it can be read by LAMMPS with the
:doc:`include <include>` command. In combination with the nocoeff
option of :doc:`write_data <write_data>` this can be used to move
the Coeffs sections from a data file into a separate file.

.. note::

   The write_coeff command is not yet fully implemented as
   some pair styles do not output their coefficient information.
   This means you will need to add/copy this information manually.

----------

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`write_restart <write_restart>`,
:doc:`write_data <write_data>`
