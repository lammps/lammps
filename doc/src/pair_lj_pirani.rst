.. index:: pair_style lj/pirani
.. index:: pair_style lj/pirani/omp

pair_style lj/pirani command
============================

Accelerator Variants: *lj/pirani/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/pirani cutoff

* lj/pirani = name of the pair style
* cutoff = global cutoff (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/pirani 10.0
   pair_coeff 1 1 4.0 7.0 6.0 3.5 0.0045

Description
"""""""""""

.. versionadded:: 12Jun2025

Pair style *lj/pirani* computes pairwise interactions from an Improved
Lennard-Jones (ILJ) potential according to :ref:`(Pirani) <Pirani>`.
The ILJ force field is adequate to model both equilibrium and
non-equilibrium properties of matter, in gaseous and condensed phases,
and at gas-surface interfaces. In particular, its use improves the
description of elementary process dynamics where the traditional
Lennard-Jones (LJ) formulation is usually applied.


.. math::

   x = r/R_m   \\
   n_x = \alpha*x^2 + \beta   \\
   \gamma \equiv m  \\

  V(x) = \varepsilon \cdot \left( \frac{\gamma}{ n_x - \gamma}  \left(\frac{1}{x} \right)^{n_x}
          -  \frac{n_x}{n_x - \gamma}  \left(\frac{1}{x} \right)^{\gamma} \right) \qquad r < r_c

:math:`r_c` is the cutoff.


An additional parameter, :math:`\alpha`, has been introduced in order to
be able to recover the traditional Lennard-Jones 12-6 with a specific
choice of parameters. With :math:`R_m \equiv r_0 = \sigma \cdot 2^{1 /
6}`, :math:`\alpha = 0`, :math:`\beta = 12` and :math:`\gamma = 6` it is
straightforward to prove that LJ 12-6 is obtained. Also, it can be
verified that using :math:`\alpha= 4`, :math:`\beta= 8` and
:math:`\gamma = 6`, at the equilibrium distance, the first and second
derivatives of ILJ match those of LJ 12-6. The parameter :math:`R_m`
corresponds to the equilibrium distance and :math:`\epsilon` to the well
depth.


This potential provides some advantages with respect to the standard LJ
potential, as explained in :ref:`(Pirani) <Pirani>`: it provides a more
realistic description of the long range behavior and an attenuation of
the hardness of the repulsive wall.

This force field can be used for neutral-neutral (:math:`\gamma = 6`),
ion-neutral (:math:`\gamma = 4`) or ion-ion systems (:math:`\gamma =
1`).  Notice that this implementation does not include explicit
electrostatic interactions.  If these are desired, this pair style
should be used along with a Coulomb pair style like
:doc:`pair styles coul/cut or coul/long <pair_coul>` by using
:doc:`pair style hybrid/overlay <pair_hybrid>` and a suitable
:doc:`kspace style <kspace_style>`, if needed.

As discussed in :ref:`(Pirani) <Pirani>`, analysis of a variety of
systems showed that :math:`\alpha= 4` generally works very well.  In
some special cases (e.g. those involving very small multiple charged
ions) this factor may take a slightly different value. The parameter
:math:`\beta` codifies the hardness (polarizability) of the interacting
partners, and for neutral-neutral systems it usually ranges from 6
to 11.  Moreover, the modulation of :math:`\beta` can model additional
interaction effects, such as charge transfer in the perturbative limit,
and can mitigate the effect of some uncertainty in the data used to
build up the potential function.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`\alpha` (dimensionless)
* :math:`\beta` (dimensionless)
* :math:`\gamma` (dimensionless)
* :math:`R_m` (distance units)
* :math:`\epsilon` (energy units)
* cutoff (distance units)

The last coefficient is optional. If not specified, the global cutoff is used.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.  Thus, coefficients for all I,J
pairs must be specified explicitly.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table options are not relevant for
this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style supports the use of the *inner*, *middle*, and
*outer* keywords of the :doc:`run_style respa <run_style>` command,
meaning the pairwise forces can be partitioned by distance at different
levels of the rRESPA hierarchy. See the :doc:`run_style <run_style>`
command for details.


----------

Restrictions
""""""""""""

This pair style is only enabled if LAMMPS was built with the EXTRA-PAIR
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`
* :doc:`pair_style lj/cut <pair_lj>`

Default
"""""""

none

--------------

.. _Pirani:

**(Pirani)** F. Pirani, S. Brizi, L. Roncaratti, P. Casavecchia, D. Cappelletti and F. Vecchiocattivi,
Phys. Chem. Chem. Phys., 2008, 10, 5489-5503.
