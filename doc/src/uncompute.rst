.. index:: uncompute

uncompute command
=================

Syntax
""""""


.. parsed-literal::

   uncompute compute-ID

* compute-ID = ID of a previously defined compute

Examples
""""""""


.. parsed-literal::

   uncompute 2
   uncompute lower-boundary

Description
"""""""""""

Delete a compute that was previously defined with a :doc:`compute <compute>`
command.  This also wipes out any additional changes made to the compute
via the :doc:`compute\_modify <compute_modify>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
