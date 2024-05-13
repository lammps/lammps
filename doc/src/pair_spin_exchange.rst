.. index:: pair_style spin/exchange
.. index:: pair_style spin/exchange/biquadratic

pair_style spin/exchange command
================================

pair_style spin/exchange/biquadratic command
============================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style spin/exchange cutoff
   pair_style spin/exchange/biquadratic cutoff

* cutoff = global cutoff pair (distance in metal units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style spin/exchange 4.0
   pair_coeff * * exchange 4.0 0.0446928 0.003496 1.4885
   pair_coeff 1 2 exchange 6.0 -0.01575 0.0 1.965 offset yes

   pair_style spin/exchange/biquadratic 4.0
   pair_coeff * * biquadratic 4.0 0.05 0.03 1.48 0.05 0.03 1.48 offset no
   pair_coeff 1 2 biquadratic 6.0 -0.01 0.0 1.9 0.0 0.1 19

Description
"""""""""""

Style *spin/exchange* computes the exchange interaction between
pairs of magnetic spins:

.. math::

   H_{ex} = -\sum_{i,j}^N J_{ij} (r_{ij}) \,\vec{s}_i \cdot \vec{s}_j

where :math:`\vec{s}_i` and :math:`\vec{s}_j` are two unit vectors representing
the magnetic spins of two particles (usually atoms), and
:math:`r_{ij} = \vert \vec{r}_i - \vec{r}_j \vert` is the inter-atomic distance
between those two particles. The summation is over pairs of nearest neighbors.
:math:`J(r_{ij})` is a function defining the intensity and the sign of the
exchange interaction for different neighboring shells.

Style *spin/exchange/biquadratic* computes a biquadratic exchange interaction
between pairs of magnetic spins:

.. math::

   H_{bi} = -\sum_{i, j}^{N} {J}_{ij} \left(r_{ij} \right)\,
                      \vec{s}_{i}\cdot \vec{s}_{j}
                      -\sum_{i, j}^{N} {K}_{ij} \left(r_{ij} \right)\,
                      \left(\vec{s}_{i}\cdot
                      \vec{s}_{j}\right)^2

where :math:`\vec{s}_i`,  :math:`\vec{s}_j`,  :math:`r_{ij}` and
:math:`J(r_{ij})` have the same definitions as above, and :math:`K(r_{ij})` is
a second function, defining the intensity and the sign of the biquadratic term.

The interatomic dependence of :math:`J(r_{ij})` and :math:`K(r_{ij})` in both
interactions above is defined by the following function:

.. math::

    {f}\left( r_{ij} \right) = 4 a \left( \frac{r_{ij}}{d}  \right)^2
    \left( 1 - b \left( \frac{r_{ij}}{d}  \right)^2 \right)
    e^{-\left( \frac{r_{ij}}{d} \right)^2 }\Theta (R_c - r_{ij})

where :math:`a`, :math:`b` and :math:`d` are the three constant coefficients
defined in the associated "pair_coeff" command, and :math:`R_c` is the radius
cutoff associated to the pair interaction (see below for more explanations).

The coefficients :math:`a`, :math:`b`, and :math:`d` need to be fitted so that
the function above matches with the value of the exchange interaction for the
:math:`N` neighbor shells taken into account.
Examples and more explanations about this function and its parameterization
are reported in :ref:`(Tranchida) <Tranchida3>`.

When a *spin/exchange/biquadratic* pair style is defined, six coefficients
(three for :math:`J(r_{ij})`, and three for :math:`K(r_{ij})`) have to be
fitted.

From this exchange interaction, each spin :math:`i` will be submitted
to a magnetic torque :math:`\vec{\omega}_{i}`, and its associated atom can be
submitted to a force :math:`\vec{F}_{i}` for spin-lattice calculations (see
:doc:`fix nve/spin <fix_nve_spin>`), such as:

.. math::

   \vec{\omega}_{i} = \frac{1}{\hbar} \sum_{j}^{Neighb} {J}
   \left(r_{ij} \right)\,\vec{s}_{j}
   ~~{\rm and}~~
   \vec{F}_{i} = \sum_{j}^{Neighb} \frac{\partial {J} \left(r_{ij} \right)}{
   \partial r_{ij}} \left( \vec{s}_{i}\cdot \vec{s}_{j} \right) \vec{e}_{ij}

with :math:`\hbar` the Planck constant (in metal units), and :math:`\vec{e}_{ij}
= \frac{\vec{r}_i - \vec{r}_j}{\vert \vec{r}_i-\vec{r}_j \vert}` the unit
vector between sites :math:`i` and :math:`j`.
Equivalent forces and magnetic torques are generated for the biquadratic term
when a *spin/exchange/biquadratic* pair style is defined.

More details about the derivation of these torques/forces are reported in
:ref:`(Tranchida) <Tranchida3>`.

For the *spin/exchange* and *spin/exchange/biquadratic* pair styles, the
following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, and set in the following order:

* :math:`R_c` (distance units)
* :math:`a`  (energy units)
* :math:`b`  (adim parameter)
* :math:`d`  (distance units)

for the *spin/exchange* pair style, and:

* :math:`R_c` (distance units)
* :math:`a_j`  (energy units)
* :math:`b_j`  (adim parameter)
* :math:`d_j`  (distance units)
* :math:`a_k`  (energy units)
* :math:`b_k`  (adim parameter)
* :math:`d_k`  (distance units)

for the *spin/exchange/biquadratic* pair style.

Note that :math:`R_c` is the radius cutoff of the considered exchange
interaction, and :math:`a`, :math:`b` and :math:`d` are the three coefficients
performing the parameterization of the function :math:`J(r_{ij})` defined
above (in the *biquadratic* style, :math:`a_j`, :math:`b_j`, :math:`d_j` and
:math:`a_k`, :math:`b_k`, :math:`d_k` are the coefficients of :math:`J(r_{ij})`
and :math:`K(r_{ij})` respectively).


None of those coefficients is optional. If not specified, the
*spin/exchange* pair style cannot be used.

----------

**Offsetting magnetic forces and energies**\ :

For spin-lattice simulation, it can be useful to offset the
mechanical forces and energies generated by the exchange
interaction.
The *offset* keyword allows to apply this offset.
By setting *offset* to *yes*, the energy definitions above are
replaced by:

.. math::

   H_{ex} = -\sum_{i,j}^N J_{ij} (r_{ij}) \,[ \vec{s}_i \cdot \vec{s}_j-1 ]

for the *spin/exchange* pair style, and:

.. math::

   H_{bi} = -\sum_{i, j}^{N} {J}_{ij} \left(r_{ij} \right)\,
                      [ \vec{s}_{i}\cdot \vec{s}_{j} -1 ]
                      -\sum_{i, j}^{N} {K}_{ij} \left(r_{ij} \right)\,
                      [ \left(\vec{s}_{i}\cdot
                      \vec{s}_{j}\right)^2 -1]

for the *spin/exchange/biquadratic* pair style.

Note that this offset only affects the calculation of the energy
and mechanical forces. It does not modify the calculation of the
precession vectors (and thus does no impact the purely magnetic
properties).
This ensures that when all spins are aligned, the magnetic energy
and the associated mechanical forces (and thus the pressure
generated by the magnetic potential) are null.

.. note::
  This offset term can be very important when calculations such as
  equations of state (energy vs volume, or energy vs pressure) are
  being performed. Indeed, setting the *offset* term ensures that
  at the ground state of the crystal and at the equilibrium magnetic
  configuration (typically ferromagnetic), the pressure is null,
  as expected.
  Otherwise, magnetic forces could generate a residual pressure.

When the *offset* option is set to *no*, no offset is applied
(also corresponding to the default option).

----------

Restrictions
""""""""""""

All the *pair/spin* styles are part of the SPIN package.  These styles
are only enabled if LAMMPS was built with this package, and if the
atom_style "spin" was declared.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`atom_style spin <atom_style>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`pair_eam <pair_eam>`,

Default
"""""""


The default *offset* keyword value is *no*.

----------

.. _Tranchida3:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
