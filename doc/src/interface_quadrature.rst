.. index:: interface_quadrature

interface_quadrature command
=============================

Syntax
""""""

.. parsed-literal::

   interface_quadrature flag

* flag = *on* or *off*

Examples
""""""""

.. code-block:: LAMMPS

   interface_quadrature off

Description
"""""""""""

Forces the force calculation process for a CAC model to use as many quadrature points as
the force cutoff range requires at the numerical interfaces between smaller elements (or atoms)
and larger elements; this can improve model stability. "Smaller" in this context refers
to elements on the other side of a numerical interface having at least 50% more undeformed volume
than elements on the former side. When this setting is on, elements
at numerical interfaces will ignore the use of the *one* keyword that can be optionally set
for a CAC :doc:`pair_style <pair_style>`.

Restrictions
""""""""""""
 This command requires a cac atom style and should be used after
reading a data or restart file.

**Default:** on
