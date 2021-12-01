.. index:: fix lb/pc

fix lb/pc command
=================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID lb/pc

* ID, group-ID are documented in the :doc:`fix <fix>` command
* lb/pc = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all lb/pc

Description
"""""""""""

This fix was part of the old LATBOLTZ package and is now defunct.  LAMMPS standard integrator :doc:`fix NVE <fix_nve>` can be used in its place.
