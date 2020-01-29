.. index:: improper\_style class2

improper\_style class2 command
==============================

improper\_style class2/omp command
==================================

improper\_style class2/kk command
=================================

Syntax
""""""


.. parsed-literal::

   improper_style class2

Examples
""""""""


.. parsed-literal::

   improper_style class2
   improper_coeff 1 100.0 0
   improper_coeff \* aa 0.0 0.0 0.0 115.06 130.01 115.06

Description
"""""""""""

The *class2* improper style uses the potential

.. image:: Eqs/improper_class2.jpg
   :align: center

where Ei is the improper term and Eaa is an angle-angle term.  The 3 X
terms in Ei are an average over 3 out-of-plane angles.

The 4 atoms in an improper quadruplet (listed in the data file read by
the :doc:`read_data <read_data>` command) are ordered I,J,K,L.  X\_IJKL
refers to the angle between the plane of I,J,K and the plane of J,K,L,
and the bond JK lies in both planes.  Similarly for X\_KJLI and X\_LJIK.
Note that atom J appears in the common bonds (JI, JK, JL) of all 3 X
terms.  Thus J (the 2nd atom in the quadruplet) is the atom of
symmetry in the 3 X angles.

The subscripts on the various theta's refer to different combinations
of 3 atoms (I,J,K,L) used to form a particular angle.  E.g. Theta\_IJL
is the angle formed by atoms I,J,L with J in the middle.  Theta1,
theta2, theta3 are the equilibrium positions of those angles.  Again,
atom J (the 2nd atom in the quadruplet) is the atom of symmetry in the
theta angles, since it is always the center atom.

Since atom J is the atom of symmetry, normally the bonds J-I, J-K, J-L
would exist for an improper to be defined between the 4 atoms, but
this is not required.

See :ref:`(Sun) <improper-Sun>` for a description of the COMPASS class2 force field.

Coefficients for the Ei and Eaa formulas must be defined for each
improper type via the :doc:`improper_coeff <improper_coeff>` command as
in the example above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands.

These are the 2 coefficients for the Ei formula:

* K (energy/radian\^2)
* X0 (degrees)

X0 is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of K are in energy/radian\^2.

For the Eaa formula, each line in a
:doc:`improper_coeff <improper_coeff>` command in the input script lists
7 coefficients, the first of which is "aa" to indicate they are
AngleAngle coefficients.  In a data file, these coefficients should be
listed under a "AngleAngle Coeffs" heading and you must leave out the
"aa", i.e. only list 6 coefficients after the improper type.

* aa
* M1 (energy/distance)
* M2 (energy/distance)
* M3 (energy/distance)
* theta1 (degrees)
* theta2 (degrees)
* theta3 (degrees)

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


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
