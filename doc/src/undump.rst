.. index:: undump

undump command
==============

Syntax
""""""

.. code-block:: LAMMPS

   undump dump-ID

* dump-ID = ID of previously defined dump

Examples
""""""""

.. code-block:: LAMMPS

   undump mine
   undump 2

Description
"""""""""""

Turn off a previously defined dump so that it is no longer active.
This closes the file associated with the dump.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dump <dump>`

Default
"""""""

none
