.. index:: pair\_style spin/neel

pair\_style spin/neel command
=============================

Syntax
""""""


.. parsed-literal::

   pair_style spin/neel cutoff

* cutoff = global cutoff pair (distance in metal units)


Examples
""""""""


.. parsed-literal::

   pair_style spin/neel 4.0
   pair_coeff \* \* neel 4.0 0.0048 0.234 1.168 2.6905 0.705 0.652
   pair_coeff 1 2 neel 4.0 0.0048 0.234 1.168 0.0 0.0 1.0

Description
"""""""""""

Style *spin/neel* computes the Neel pair anisotropy model
between pairs of magnetic spins:

.. math source doc: src/Eqs/pair_spin_neel_interaction.tex
.. math::

   :align: center

where si and sj are two neighboring magnetic spins of two particles,
rij = ri - rj is the inter-atomic distance between the two particles,
eij = (ri - rj)/\|ri-rj\| is their normalized separation vector and g1,
q1 and q2 are three functions defining the intensity of the dipolar
and quadrupolar contributions, with:

.. math source doc: src/Eqs/pair_spin_neel_functions.tex
.. math::

   :align: center

With the functions g(rij) and q(rij) defined and fitted according to
the same Bethe-Slater function used to fit the exchange interaction:

.. math source doc: src/Eqs/pair_spin_exchange_function.tex
.. math::

   :align: center

where a, b and d are the three constant coefficients defined in the
associated "pair\_coeff" command.

The coefficients a, b, and d need to be fitted so that the function
above matches with the values of the magneto-elastic constant of the
materials at stake.

Examples and more explanations about this function and its
parameterization are reported in :ref:`(Tranchida) <Tranchida6>`. More
examples of parameterization will be provided in future work.

From this DM interaction, each spin i will be submitted to a magnetic
torque omega and its associated atom to a force F (for spin-lattice
calculations only).

More details about the derivation of these torques/forces are reported
in :ref:`(Tranchida) <Tranchida6>`.


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


.. _Tranchida6:



**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
