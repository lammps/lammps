.. index:: compute msd/nongauss

compute msd/nongauss command
============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID msd/nongauss keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* msd/nongauss = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *com*

  .. parsed-literal::

       *com* value = *yes* or *no*

Examples
""""""""

.. parsed-literal::

   compute 1 all msd/nongauss
   compute 1 upper msd/nongauss com yes

Description
"""""""""""

Define a computation that calculates the mean-squared displacement
(MSD) and non-Gaussian parameter (NGP) of the group of atoms,
including all effects due to atoms passing through periodic boundaries.

A vector of three quantities is calculated by this compute.  The first
element of the vector is the total squared dx,dy,dz displacements
drsquared = (dx\*dx + dy\*dy + dz\*dz) of atoms, and the second is the
fourth power of these displacements drfourth = (dx\*dx + dy\*dy +
dz\*dz)\*(dx\*dx + dy\*dy + dz\*dz), summed and averaged over atoms in the
group.  The 3rd component is the nonGaussian diffusion parameter NGP =
3\*drfourth/(5\*drsquared\*drsquared), i.e.

.. math::

 NGP(t) = 3<(r(t)-r(0))^4>/(5<(r(t)-r(0))^2>^2) - 1

The NGP is a commonly used quantity in studies of dynamical
heterogeneity.  Its minimum theoretical value (-0.4) occurs when all
atoms have the same displacement magnitude.  NGP=0 for Brownian
diffusion, while NGP > 0 when some mobile atoms move faster than
others.

If the *com* option is set to *yes* then the effect of any drift in
the center-of-mass of the group of atoms is subtracted out before the
displacement of each atom is calculated.

See the :doc:`compute msd <compute_msd>` doc page for further important
NOTEs, which also apply to this compute.

**Output info:**

This compute calculates a global vector of length 3, which can be
accessed by indices 1-3 by any command that uses global vector values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The vector values are "intensive".  The first vector value will be in
distance\^2 :doc:`units <units>`, the second is in distance\^4 units, and
the 3rd is dimensionless.

Restrictions
""""""""""""

This compute is part of the MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute msd <compute_msd>`

Default
"""""""

The option default is com = no.
