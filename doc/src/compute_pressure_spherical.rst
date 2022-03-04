.. index:: compute pressure/spherical

compute pressure/spherical command
==================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID pressure/spherical x0 y0 z0 bin_width Rmax

* ID, group-ID are documented in :doc:`compute <compute>` command
* pressure/cartesian = style name of this compute command
* x0, y0, z0 = origin of the spherical coordinate system
* bin_width = width of spherical shells
* Rmax = maximum radius of spherical shells

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all pressure/spherical 0 0 0 0.1 10

Description
"""""""""""

This computes the diagonal components of the spherical pressure tensor in spherical shells, as described in :ref:`(Ikeshoji)<Ikeshoji3>`. This can be used to calculate the local pressure tensor components of for example droplets and bubbles. This compute obeys momentum balance. This compute uses the Irving-Kirkwood contour, whici is the straight line between particle pairs. The pressure tensor is split into a kinetic contribution :math:`P^k` and a virial contribution :math:`P^v`. The sum gives the pressure :math:`P = P^k+P^v`.

Output info
"""""""""""

This compute outputs a global array with 8 columns and Rmax/bin_width. The output columns are position of the center of the local spherical shell, number density, :math:`P^k_{rr}`, :math:`P^k_{\theta\theta}`, :math:`P^k_{\phi\phi}`, :math:`P^v_{rr}`, :math:`P^v_{\theta\theta}`, and :math:`P^v_{\phi\phi}`.

This array can be output with :doc:`fix ave/time <fix_ave_time>`,

.. code-block:: LAMMPS

  compute p all pressure/spherical 0 0 0 0.1 10
  fix 2 all ave/time 100 1 100 c_p[*] file dump_p.out mode vector

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute stress/atom <compute_stress_atom>`, :doc:`compute pressure/cylinder <compute_pressure_cylinder>`, :doc:`compute pressure/cartesian <compute_pressure_cartesian>`

Default
"""""""

none

----------

.. _Ikeshoji3:
**(Ikeshoji)** Ikeshoji, Hafskjold, Furuholt, Mol Sim, 29, 101-109, (2003).
