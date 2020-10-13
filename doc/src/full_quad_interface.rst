.. index:: full_quad_interface

full_quad_interface command
=============================

Syntax
""""""

.. parsed-literal::

   full_quad_interface flag

* flag = *on* or *off*

Examples
""""""""

.. code-block:: LAMMPS

   full_quad_interface on

Description
"""""""""""

Force the force calculation process for a CAC model to use more quadrature points at numerical
interfaces between smaller elements (or atoms) and larger elements; this can improve model stability.
Currently, this command loops over all possible atom sites of elements at numerical interfaces; thus,
refrain from constructing models with sudden transitions in element length scales.

Restrictions
""""""""""""
 This command requires a cac atom style and should be used after reading a data or restart file.

**Default:** off
