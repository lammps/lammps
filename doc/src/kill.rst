.. index:: kill

kill command
============

Syntax
""""""

.. code-block:: LAMMPS

   kill rank kill_rank

* kill_rank = rank to kill

Examples
""""""""

.. code-block:: LAMMPS

   kill rank 1

Description
"""""""""""

Used for testing fault handling. Generates an artificial failure on kill_rank
by raising the SIGTERM signal. Calls std::abort() if SIGTERM is ignored.

Restrictions
""""""""""""

This fix is part of the FENIX package. It is only enabled if LAMMPS was built
with that package. See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`fix fail` :doc:`fenix`
