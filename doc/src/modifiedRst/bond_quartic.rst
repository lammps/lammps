.. index:: bond\_style quartic

bond\_style quartic command
===========================

bond\_style quartic/omp command
===============================

Syntax
""""""


.. parsed-literal::

   bond_style quartic

Examples
""""""""


.. parsed-literal::

   bond_style quartic
   bond_coeff 2 1200 -0.55 0.25 1.3 34.6878

Description
"""""""""""

The *quartic* bond style uses the potential

.. math::

  E = K (r - R_c)^ 2 (r - R_c - B_1) (r - R_c - B_2) + U_0 +
  4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - 
    \left(\frac{\sigma}{r}\right)^6 \right] + \epsilon


to define a bond that can be broken as the simulation proceeds (e.g.
due to a polymer being stretched).  The sigma and epsilon used in the
LJ portion of the formula are both set equal to 1.0 by LAMMPS.

The following coefficients must be defined for each bond type via the
:doc:`bond\_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* K (energy/distance\^4)
* B1 (distance)
* B2 (distance)
* Rc (distance)
* U0 (energy)

This potential was constructed to mimic the FENE bond potential for
coarse-grained polymer chains.  When monomers with sigma = epsilon =
1.0 are used, the following choice of parameters gives a quartic
potential that looks nearly like the FENE potential: K = 1200, B1 =
-0.55, B2 = 0.25, Rc = 1.3, and U0 = 34.6878.  Different parameters
can be specified using the :doc:`bond\_coeff <bond_coeff>` command, but
you will need to choose them carefully so they form a suitable bond
potential.

Rc is the cutoff length at which the bond potential goes smoothly to a
local maximum.  If a bond length ever becomes > Rc, LAMMPS "breaks"
the bond, which means two things.  First, the bond potential is turned
off by setting its type to 0, and is no longer computed.  Second, a
pairwise interaction between the two atoms is turned on, since they
are no longer bonded.

LAMMPS does the second task via a computational sleight-of-hand.  It
subtracts the pairwise interaction as part of the bond computation.
When the bond breaks, the subtraction stops.  For this to work, the
pairwise interaction must always be computed by the
:doc:`pair\_style <pair_style>` command, whether the bond is broken or
not.  This means that :doc:`special\_bonds <special_bonds>` must be set
to 1,1,1, as indicated as a restriction below.

Note that when bonds are dumped to a file via the :doc:`dump local <dump>` command, bonds with type 0 are not included.  The
:doc:`delete\_bonds <delete_bonds>` command can also be used to query the
status of broken bonds or permanently delete them, e.g.:


.. parsed-literal::

   delete_bonds all stats
   delete_bonds all bond 0 remove


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


This bond style can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

The *quartic* style requires that :doc:`special\_bonds <special_bonds>`
parameters be set to 1,1,1.  Three- and four-body interactions (angle,
dihedral, etc) cannot be used with *quartic* bonds.

Related commands
""""""""""""""""

:doc:`bond\_coeff <bond_coeff>`, :doc:`delete\_bonds <delete_bonds>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
