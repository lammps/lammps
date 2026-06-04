.. index:: pair_style lj/cut/tip4p/long/functor

pair_style lj/cut/tip4p/long/functor command
============================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cut/tip4p/long/functor typeO typeH typeB typeA qdist cutoff (cutoff2)

* typeO,typeH = atom types of TIP4P oxygen and hydrogen
* typeB,typeA = bond and angle types of the TIP4P molecule
* qdist = distance from O atom to the M charge site (distance units)
* cutoff = global cutoff for LJ (and Coulomb if only one cutoff) (distance units)
* cutoff2 = global cutoff for Coulomb (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/tip4p/long/functor 1 2 1 1 0.15 12.0
   pair_coeff * * 0.155 3.15
   pair_coeff 2 2 0.0 0.0

Description
"""""""""""

.. versionadded:: TBD

Style *lj/cut/tip4p/long/functor* computes the cutoff Lennard-Jones potential
together with a long-range (Ewald/PPPM) Coulomb term in which the charge of the
TIP4P water oxygen is placed on a virtual *M* site, and is numerically
equivalent to :doc:`pair_style lj/cut/tip4p/long <pair_lj_cut_tip4p>`.  It is
provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation using
the template/functor framework described in :doc:`Developer_write_pair_functor`
(the Lennard-Jones evaluator combined with a dedicated TIP4P driver that reuses
the long-range-Coulomb policy for the real-space electrostatics).

The arguments and ``pair_coeff`` coefficients are as in
:doc:`pair_style lj/cut/tip4p/long <pair_lj_cut_tip4p>`.  It must be used with a
:doc:`kspace_style <kspace_style>` that supports TIP4P (e.g. *pppm/tip4p*) and
with bond and angle styles (for the M-site geometry).

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.  It additionally requires a TIP4P-capable
:doc:`KSpace style <kspace_style>` (from the KSPACE package), atom IDs,
:doc:`newton <newton>` *on* for pairwise interactions, and a charge per atom.

The Coulomb cutoff is global.

Related commands
""""""""""""""""

:doc:`pair_style lj/cut/tip4p/long <pair_lj_cut_tip4p>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
