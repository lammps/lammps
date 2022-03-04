.. index:: compute pressure/cartesian

compute pressure/cartesian command
==================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID pressure/cartesian dim bin_width 

* ID, group-ID are documented in :doc:`compute <compute>` command
* pressure/cartesian = style name of this compute command
* one or two dim/bin_width pairs may be appended
* dim = x, y, or z
* bin_width = width of the bin 

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all pressure/cartesian x 0.1
   compute 1 all pressure/cartesian y 0.25 z 0.1

Description
"""""""""""

This computes the Cartesian pressure tensor in one or two dimensions, as described in :ref:`(Ikeshoji)<Ikeshoji2>`. This can for example be used to calculate the local pressure tensor of flat liquid-vapor interfaces. This compute obeys momentum balance, such that the pressure tensor component normal to a surface is constant. This compute uses the Irving-Kirkwood contour, which is the straight line between particle pairs. The pressure tensor is split into a kinetic contribution :math:`P^k` and a virial contribution :math:`P^v`. The sum gives the total pressure tensor :math:`P = P^k+P^v`.

Output info
"""""""""""

The output is a global array with 8 columns when only one dimension is specified and 9 columns when two dimensions are specified. There are (L1*L2)/(bin_width1*bin_width2) rows, where L1 and L2 are the size of the simulation box in the relevant dimensions. The output columns are position of the center of the local volume in the first and second dimension (when two dimensions are specified), number density, :math:`P^k_{xx}`, :math:`P^k_{yy}`, :math:`P^k_{zz}`, :math:`P^v_{xx}`, :math:`P^v_{yy}`, and :math:`P^v_{zz}`.

This array can be output with the :doc:`fix ave/time <fix_ave_time>`,

.. code-block:: LAMMPS

  compute p all pressure/cartesian x 0.1
  fix 2 all ave/time 100 1 100 c_p[*] file dump_p.out mode vector

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute stress/atom <compute_stress_atom>`, :doc:`compute pressure/cylinder <compute_pressure_cylinder>`, :doc:`compute pressure/spherical <compute_pressure_spherical>`

Default
"""""""

none

----------

.. _Ikeshoji2:

**(Ikeshoji)** Ikeshoji, Hafskjold, Furuholt, Mol Sim, 29, 101-109, (2003).
