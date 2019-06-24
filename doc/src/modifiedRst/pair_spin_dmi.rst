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

.. math source doc: src/Eqs/pair_spin_dmi_interaction.tex
.. math::

   :align: center

where si and sj are two neighboring magnetic spins of two particles,
eij = (ri - rj)/\|ri-rj\| is the unit vector between sites i and j,
and D is the DM vector defining the intensity (in eV) and the direction
of the interaction.

In :ref:`(Rohart) <Rohart>`, D is defined as the direction normal to the film oriented
from the high spin-orbit layer to the magnetic ultra-thin film.

The application of a spin-lattice Poisson bracket to this energy (as described
in :ref:`(Tranchida) <Tranchida5>`) allows to derive a magnetic torque omega, and a
mechanical force F (for spin-lattice calculations only) for each magnetic
particle i:

.. math source doc: src/Eqs/pair_spin_dmi_forces.tex
.. math::

   :align: center

More details about the derivation of these torques/forces are reported in
:ref:`(Tranchida) <Tranchida5>`.

For the *spin/dmi* pair style, the following coefficients must be defined for
each pair of atoms types via the :doc:`pair\_coeff <pair_coeff>` command as in
the examples above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>` commands, and
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

:doc:`atom\_style spin <atom_style>`, :doc:`pair\_coeff <pair_coeff>`,
:doc:`pair\_eam <pair_eam>`,

**Default:** none


----------


.. _Rohart:



.. _Tranchida5:

**(Rohart)** Rohart and Thiaville,
Physical Review B, 88(18), 184422. (2013).


**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
