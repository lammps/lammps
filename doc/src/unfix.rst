.. index:: unfix

unfix command
=============

Syntax
""""""


.. parsed-literal::

   unfix fix-ID

* fix-ID = ID of a previously defined fix

Examples
""""""""


.. parsed-literal::

   unfix 2
   unfix lower-boundary

Description
"""""""""""

Delete a fix that was previously defined with a :doc:`fix <fix>`
command.  This also wipes out any additional changes made to the fix
via the :doc:`fix\_modify <fix_modify>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix <fix>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
