.. index:: improper\_style distharm

improper\_style distharm command
================================

Syntax
""""""

improper\_style distharm

Examples
""""""""


.. parsed-literal::

   improper_style distharm
   improper_coeff 1 25.0 0.5

Description
"""""""""""

The *distharm* improper style uses the potential

.. math source doc: src/Eqs/improper_distharm.tex
.. math::

   E = K (d - d_0)^2


where d is the oriented distance between the central atom and the plane formed
by the other three atoms.  If the 4 atoms in an improper quadruplet
(listed in the data file read by the :doc:`read\_data <read_data>`
command) are ordered I,J,K,L then the L-atom is assumed to be the
central atom. Note that this is different from the convention used
in the improper\_style distance. The distance d is oriented and can take
on negative values. This may lead to unwanted behavior if d0 is not equal to zero.

The following coefficients must be defined for each improper type via
the improper\_coeff command as in the example above, or in the data
file or restart files read by the read\_data or read\_restart commands:

* K (energy/distance\^2)
* d0 (distance)


----------


Restrictions
""""""""""""


This improper style can only be used if LAMMPS was built with the
USER-YAFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`improper\_coeff <improper_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
