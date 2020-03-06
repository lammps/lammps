.. index:: improper_style class2

improper_style class2 command
=============================

improper_style class2/omp command
=================================

improper_style class2/kk command
================================

Syntax
""""""


.. code-block:: LAMMPS

   improper_style class2

Examples
""""""""


.. code-block:: LAMMPS

   improper_style class2
   improper_coeff 1 100.0 0
   improper_coeff * aa 0.0 0.0 0.0 115.06 130.01 115.06

Description
"""""""""""

The *class2* improper style uses the potential

.. math::

   E      = & E_i + E_{aa} \\
   E_i    = & K [ \frac{\chi_{ijkl} + \chi_{kjli} + \chi_{ljik}}{3} - \chi_0 ]^2 \\
   E_{aa} = & M_1 (\theta_{ijk} - \theta_1) (\theta_{kjl} - \theta_3) + \\
            & M_2 (\theta_{ijk} - \theta_1) (\theta_{ijl} - \theta_2) + \\
            & M_3 (\theta_{ijl} - \theta_2) (\theta_{kjl} - \theta_3)


where :math:`E_i` is the improper term and :math:`E_{aa}` is an
angle-angle term.  The 3 :math:`\chi` terms in :math:`E_i` are an
average over 3 out-of-plane angles.

The 4 atoms in an improper quadruplet (listed in the data file read by
the :doc:`read_data <read_data>` command) are ordered I,J,K,L.
:math:`\chi_{ijkl}` refers to the angle between the plane of I,J,K and
the plane of J,K,L, and the bond JK lies in both planes.  Similarly for
:math:`\chi_{kjli}` and :math:`\chi_{ljik}`.
Note that atom J appears in the common bonds (JI, JK, JL) of all 3 X
terms.  Thus J (the 2nd atom in the quadruplet) is the atom of
symmetry in the 3 :math:`\chi` angles.

The subscripts on the various :math:`\theta`\ s refer to different
combinations of 3 atoms (I,J,K,L) used to form a particular angle.
E.g. :math:`\theta_{ijl}` is the angle formed by atoms I,J,L with J
in the middle.  :math:`\theta_1`, :math:`\theta_2`, :math:`\theta_3`
are the equilibrium positions of those angles.  Again,
atom J (the 2nd atom in the quadruplet) is the atom of symmetry in the
theta angles, since it is always the center atom.

Since atom J is the atom of symmetry, normally the bonds J-I, J-K, J-L
would exist for an improper to be defined between the 4 atoms, but
this is not required.

See :ref:`(Sun) <improper-Sun>` for a description of the COMPASS class2 force field.

Coefficients for the :math:`E_i` and :math:`E_{aa}` formulas must be
defined for each
improper type via the :doc:`improper_coeff <improper_coeff>` command as
in the example above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands.

These are the 2 coefficients for the :math:`E_i` formula:

* :math:`K` (energy/radian\^2)
* :math:`\chi_0` (degrees)

:math:`\chi_0` is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of K are in energy/radian\^2.

For the :math:`E_{aa}` formula, each line in a
:doc:`improper_coeff <improper_coeff>` command in the input script lists
7 coefficients, the first of which is *aa* to indicate they are
AngleAngle coefficients.  In a data file, these coefficients should be
listed under a *AngleAngle Coeffs* heading and you must leave out the
*aa*, i.e. only list 6 coefficients after the improper type.

* *aa*
* :math:`M_1` (energy/distance)
* :math:`M_2` (energy/distance)
* :math:`M_3` (energy/distance)
* :math:`\theta_1` (degrees)
* :math:`\theta_2` (degrees)
* :math:`\theta_3` (degrees)

The theta values are specified in degrees, but LAMMPS converts them to
radians internally; hence the units of M are in energy/radian\^2.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This improper style can only be used if LAMMPS was built with the
CLASS2 package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

**Default:** none


----------


.. _improper-Sun:



**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).
