.. index:: uncompute

uncompute command
=================

Syntax
""""""

.. code-block:: LAMMPS

   uncompute compute-ID

* compute-ID = ID of a previously defined compute

Examples
""""""""

.. code-block:: LAMMPS

   uncompute 2
   uncompute lower-boundary

Description
"""""""""""

Delete a compute that was previously defined with a :doc:`compute <compute>`
command.  This also wipes out any additional changes made to the compute
via the :doc:`compute_modify <compute_modify>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`

Default
"""""""

none
