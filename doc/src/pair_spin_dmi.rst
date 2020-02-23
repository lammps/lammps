.. index:: pair\_style spin/dmi

pair\_style spin/dmi command
============================

Syntax
""""""


.. parsed-literal::

   pair_style spin/dmi cutoff

* cutoff = global cutoff pair (distance in metal units)


Examples
""""""""


.. parsed-literal::

   pair_style spin/dmi 4.0
   pair_coeff \* \* dmi 2.6 0.001 1.0 0.0 0.0
   pair_coeff 1 2 dmi 4.0 0.00109 0.0 0.0 1.0

Description
"""""""""""

Style *spin/dmi* computes the Dzyaloshinskii-Moriya (DM) interaction
between pairs of magnetic spins.
According to the expression reported in :ref:`(Rohart) <Rohart>`, one has
the following DM energy:

.. math::

    \mathbf{H}_{dm} = \sum_{{ i,j}=1,i\neq j}^{N} 
    \left( \vec{e}_{ij} \times \vec{D} \right)
    \cdot\left(\vec{s}_{i}\times \vec{s}_{j}\right), 

where :math:`\vec{s}_i` and :math:`\vec{s}_j` are two neighboring magnetic spins of
two particles, :math:`\vec{e}_ij = \frac{r_i - r_j}{\left| r_i - r_j \right|}`
is the unit vector between sites *i* and *j*, and :math:`\vec{D}` is the
DM vector defining the intensity (in eV) and the direction of the
interaction.

In :ref:`(Rohart) <Rohart>`, :math:`\vec{D}` is defined as the direction normal to the film oriented
from the high spin-orbit layer to the magnetic ultra-thin film.

The application of a spin-lattice Poisson bracket to this energy (as described
in :ref:`(Tranchida) <Tranchida5>`) allows to derive a magnetic torque omega, and a
mechanical force F (for spin-lattice calculations only) for each magnetic
particle i:

.. math::

    \vec{\omega}_i = -\frac{1}{\hbar} \sum_{j}^{Neighb} \vec{s}_{j}\times \left(\vec{e}_{ij}\times \vec{D} \right) 
    ~~{\rm and}~~
    \vec{F}_i = -\sum_{j}^{Neighb} \frac{1}{r_{ij}} \vec{D} \times \left( \vec{s}_{i}\times \vec{s}_{j} \right) 

More details about the derivation of these torques/forces are reported in
:ref:`(Tranchida) <Tranchida5>`.

For the *spin/dmi* pair style, the following coefficients must be defined for
each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in
the examples above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands, and
set in the following order:

* rc (distance units)
* \|D\| (energy units)
* Dx, Dy, Dz  (direction of D)

Note that rc is the radius cutoff of the considered DM interaction, \|D\| is
the norm of the DM vector (in eV), and Dx, Dy and Dz define its direction.

None of those coefficients is optional.  If not specified, the *spin/dmi*
pair style cannot be used.


----------


Restrictions
""""""""""""


All the *pair/spin* styles are part of the SPIN package.  These styles
are only enabled if LAMMPS was built with this package, and if the
atom\_style "spin" was declared.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`atom_style spin <atom_style>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`pair_eam <pair_eam>`,

**Default:** none


----------

.. _Rohart:

.. _Tranchida5:

**(Rohart)** Rohart and Thiaville,
Physical Review B, 88(18), 184422. (2013).


**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
