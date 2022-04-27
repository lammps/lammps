.. index:: pair_style lebedeva/z

pair_style lebedeva/z command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style [hybrid/overlay ...] lebedeva/z cutoff

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay lebedeva/z 20.0
   pair_coeff * * none
   pair_coeff 1 2 lebedeva/z  CC.Lebedeva   C C

   pair_style hybrid/overlay rebo lebedeva/z 14.0
   pair_coeff * * rebo        CH.rebo       C C
   pair_coeff 1 2 lebedeva/z  CC.Lebedeva   C C

Description
"""""""""""

The *lebedeva/z* pair style computes the Lebedeva interaction potential
as described in :ref:`(Lebedeva1) <Leb01>` and :ref:`(Lebedeva2)
<Leb02>`.  An important simplification is made, which is to take all
normals along the z-axis.

The Lebedeva potential is intended for the description of the interlayer
interaction between graphene layers.  To perform a realistic simulation,
this potential must be used in combination with an intralayer potential
such as :doc:`AIREBO <pair_airebo>` or :doc:`Tersoff <pair_tersoff>`
facilitated by using pair style :doc:`hybrid/overlay <pair_hybrid>`.  To
keep the intralayer properties unaffected, the interlayer interaction
within the same layers should be avoided.  This can be achieved by
assigning different atom types to atoms of different layers (e.g. 1 and
2 in the examples above).

Other interactions can be set to zero using pair_style *none*\ .


.. math::

   E       = & \frac{1}{2} \sum_i \sum_{j \neq i} V_{ij}\\
   V_{ij}  = & B e^{-\alpha(r_{ij} - z_0)} \\
             & + C(1 + D_1\rho^2_{ij} + D_2\rho^4_{ij}) e^{-\lambda_1\rho^2_{ij}} e^{-\lambda_2 (z^2_{ij} - z^2_0)} \\
             & - A \left(\frac{z_0}{r_ij}\right)^6 + A \left( \frac{z_0}{r_c} \right)^6 \\
   \rho^2_{ij} = & x^2_{ij} + y^2_{ij} \qquad (\mathbf{n_i} \equiv \mathbf{\hat{z}})

It is important to have a sufficiently large cutoff to ensure smooth forces.
Energies are shifted so that they go continuously to zero at the cutoff assuming
that the exponential part of :math:`V_{ij}` (first term) decays sufficiently fast.
This shift is achieved by the last term in the equation for :math:`V_{ij}` above.

The provided parameter file (CC.Lebedeva) contains two sets of parameters.

- The first set (element name "C") is suitable for normal conditions and
  is taken from :ref:`(Popov1) <Popov>`
- The second set (element name "C1") is suitable for high-pressure
  conditions and is taken from :ref:`(Koziol1) <Koziol>`

Both sets contain an additional parameter, *S*, that can be used to
facilitate scaling of energies and is set to 1.0 by default.

Restrictions
""""""""""""

This pair style is part of the INTERLAYER package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`,
:doc:`pair_style none <pair_none>`,
:doc:`pair_style hybrid/overlay <pair_hybrid>`,
:doc:`pair_style drip <pair_drip>`,
:doc:`pair_style ilp/graphene/hbd <pair_ilp_graphene_hbn>`,
:doc:`pair_style kolmogorov/crespi/z <pair_kolmogorov_crespi_z>`,
:doc:`pair_style kolmogorov/crespi/full <pair_kolmogorov_crespi_full>`.

Default
"""""""

none

----------

.. _Leb01:

**(Lebedeva1)** I. V. Lebedeva, A. A. Knizhnik, A. M. Popov, Y. E. Lozovik, B. V. Potapkin, Phys. Rev. B, 84, 245437 (2011)

.. _Leb02:

**(Lebedeva2)** I. V. Lebedeva, A. A. Knizhnik, A. M. Popov, Y. E. Lozovik, B. V. Potapkin, Physica E: 44, 949-954 (2012)

.. _Popov:

**(Popov1)** A.M. Popov, I. V. Lebedeva, A. A. Knizhnik, Y. E. Lozovik and B. V. Potapkin, Chem. Phys. Lett. 536, 82-86 (2012).

.. _Koziol:

**(Koziol1)** Z. Koziol, G. Gawlik and J. Jagielski, Chinese Phys. B 28, 096101 (2019).
