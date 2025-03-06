.. index:: unfix

unfix command
=============

Syntax
""""""

.. code-block:: LAMMPS

   unfix fix-ID

* fix-ID = ID of a previously defined fix

Examples
""""""""

.. code-block:: LAMMPS

   unfix 2
   unfix lower-boundary

Description
"""""""""""

Delete a fix that was previously defined with a :doc:`fix <fix>`
command.  This also wipes out any additional changes made to the fix
via the :doc:`fix_modify <fix_modify>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix <fix>`

Default
"""""""

none
