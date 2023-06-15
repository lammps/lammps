.. index:: compute msd/nongauss

compute msd/nongauss command
============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID msd/nongauss keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* msd/nongauss = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *com*

  .. parsed-literal::

       *com* value = *yes* or *no*

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all msd/nongauss
   compute 1 upper msd/nongauss com yes

Description
"""""""""""

Define a computation that calculates the mean-squared displacement
(MSD) and non-Gaussian parameter (NGP) of the group of atoms,
including all effects due to atoms passing through periodic boundaries.

A vector of three quantities is calculated by this compute.  The first
element of the vector is the total squared displacement,
:math:`dr^2 = dx^2 + dy^2 + dz^2`, of the atoms, and the second is the
fourth power of these displacements, :math:`dr^4 = (dx^2 + dy^2 + dz^2)^2`,
summed and averaged over atoms in the group.  The third component is the
non-Gaussian diffusion parameter NGP,

.. math::

   \text{NGP}(t) = \frac{3\left\langle(r(t)-r(0))^4\right\rangle}
                        {5\left\langle(r(t)-r(0))^2\right\rangle^2} - 1.

The NGP is a commonly used quantity in studies of dynamical
heterogeneity.  Its minimum theoretical value :math:`(-0.4)` occurs when all
atoms have the same displacement magnitude.  :math:`\text{NGP}=0` for Brownian
diffusion, while :math:`\text{NGP} > 0` when some mobile atoms move faster than
others.

If the *com* option is set to *yes* then the effect of any drift in
the center-of-mass of the group of atoms is subtracted out before the
displacement of each atom is calculated.

See the :doc:`compute msd <compute_msd>` page for further important
NOTEs, which also apply to this compute.

Output info
"""""""""""

This compute calculates a global vector of length 3, which can be
accessed by indices 1--3 by any command that uses global vector values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The vector values are "intensive".  The first vector value will be in
distance\ :math:`^2` :doc:`units <units>`, the second is in
distance\ :math:`^4` units, and the third is dimensionless.

Restrictions
""""""""""""

Compute *msd/nongauss* cannot be used with a dynamic group.

This compute is part of the EXTRA-COMPUTE package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute msd <compute_msd>`

Default
"""""""

The option default is com = no.
