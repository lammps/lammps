.. index:: angle\_style cosine/buck6d

angle\_style cosine/buck6d command
==================================

Syntax
""""""


.. parsed-literal::

   angle_style cosine/buck6d

Examples
""""""""


.. parsed-literal::

   angle_style cosine/buck6d
   angle_coeff 1  cosine/buck6d  1.978350  4  180.000000

Description
"""""""""""

The *cosine/buck6d* angle style uses the potential

.. image:: Eqs/angle_cosine_buck6d.jpg
   :align: center

where K is the energy constant, n is the periodic multiplicity and
Theta0 is the equilibrium angle.

The coefficients must be defined for each angle type via the
:doc:`angle\_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands in the following order:

* K (energy)
* n
* Theta0 (degrees)

Theta0 is specified in degrees, but LAMMPS converts it to radians
internally.

Additional to the cosine term the *cosine/buck6d* angle style computes
the short range (vdW) interaction belonging to the
:doc:`pair\_buck6d <pair_buck6d_coul_gauss>` between the end atoms of the
angle.  For this reason this angle style only works in combination
with the :doc:`pair\_buck6d <pair_buck6d_coul_gauss>` styles and needs
the :doc:`special\_bonds <special_bonds>` 1-3 interactions to be weighted
0.0 to prevent double counting.


----------


Restrictions
""""""""""""


*cosine/buck6d* can only be used in combination with the
:doc:`pair\_buck6d <pair_buck6d_coul_gauss>` style and with a
:doc:`special\_bonds <special_bonds>` 0.0 weighting of 1-3 interactions.

This angle style can only be used if LAMMPS was built with the
USER-MOFFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle\_coeff <angle_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
