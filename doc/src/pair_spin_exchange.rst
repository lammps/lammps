.. index:: pair_style spin/exchange

pair_style spin/exchange command
================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style spin/exchange cutoff

* cutoff = global cutoff pair (distance in metal units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style spin/exchange 4.0
   pair_coeff * * exchange 4.0 0.0446928 0.003496 1.4885
   pair_coeff 1 2 exchange 6.0 -0.01575 0.0 1.965

Description
"""""""""""

Style *spin/exchange* computes the exchange interaction between
pairs of magnetic spins:

.. math::

   H_{ex} = -\sum_{i,j}^N J_{ij} (r_{ij}) \,\vec{s}_i \cdot \vec{s}_j

where :math:`\vec{s}_i` and :math:`\vec{s}_j` are two neighboring magnetic spins of two particles,
:math:`r_{ij} = \vert \vec{r}_i - \vec{r}_j \vert` is the inter-atomic distance between the two
particles. The summation is over pairs of nearest neighbors.
:math:`J(r_{ij})` is a function defining the intensity and the sign of the exchange
interaction for different neighboring shells. This function is defined as:

.. math::

    {J}\left( r_{ij} \right) = 4 a \left( \frac{r_{ij}}{d}  \right)^2 \left( 1 - b \left( \frac{r_{ij}}{d}  \right)^2 \right) e^{-\left( \frac{r_{ij}}{d} \right)^2 }\Theta (R_c - r_{ij})

where :math:`a`, :math:`b` and :math:`d` are the three constant coefficients defined in the associated
"pair_coeff" command, and :math:`R_c` is the radius cutoff associated to
the pair interaction (see below for more explanations).

The coefficients :math:`a`, :math:`b`, and :math:`d` need to be fitted so that the function above matches with
the value of the exchange interaction for the :math:`N` neighbor shells taken into account.
Examples and more explanations about this function and its parameterization are reported
in :ref:`(Tranchida) <Tranchida3>`.

From this exchange interaction, each spin :math:`i` will be submitted
to a magnetic torque :math:`\vec{\omega}`, and its associated atom can be submitted to a
force :math:`\vec{F}` for spin-lattice calculations (see :doc:`fix nve/spin <fix_nve_spin>`),
such as:

.. math::

   \vec{\omega}_{i} = \frac{1}{\hbar} \sum_{j}^{Neighb} {J}
   \left(r_{ij} \right)\,\vec{s}_{j}
   ~~{\rm and}~~
   \vec{F}_{i} = \sum_{j}^{Neighb} \frac{\partial {J} \left(r_{ij} \right)}{ \partial r_{ij}} \left( \vec{s}_{i}\cdot \vec{s}_{j} \right) \vec{e}_{ij}

with :math:`\hbar` the Planck constant (in metal units), and :math:`\vec{e}_{ij} = \frac{\vec{r}_i - \vec{r}_j}{\vert \vec{r}_i-\vec{r}_j \vert}` the unit
vector between sites :math:`i` and :math:`j`.

More details about the derivation of these torques/forces are reported in
:ref:`(Tranchida) <Tranchida3>`.

For the *spin/exchange* pair style, the following coefficients must be defined
for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in
the examples above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands, and
set in the following order:

* :math:`R_c` (distance units)
* :math:`a`  (energy units)
* :math:`b`  (adim parameter)
* :math:`d`  (distance units)

Note that :math:`R_c` is the radius cutoff of the considered exchange interaction,
and :math:`a`, :math:`b` and :math:`d` are the three coefficients performing the parameterization
of the function :math:`J(r_{ij})` defined above.

None of those coefficients is optional. If not specified, the
*spin/exchange* pair style cannot be used.

----------

Restrictions
""""""""""""

All the *pair/spin* styles are part of the SPIN package.  These styles
are only enabled if LAMMPS was built with this package, and if the
atom_style "spin" was declared.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`atom_style spin <atom_style>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`pair_eam <pair_eam>`,

**Default:**

none

----------

.. _Tranchida3:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
