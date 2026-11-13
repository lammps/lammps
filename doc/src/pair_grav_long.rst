.. index:: pair_style grav/long

pair_style grav/long command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style grav/long cutoff keyword value

* cutoff = global cutoff for the real-space part of the interaction (distance units)
* zero or one keyword/value pair may be appended
* keyword = *gconst*

  .. parsed-literal::

       *gconst* value = G
         G = gravitational constant (energy*distance/mass\^2 units)

Examples
""""""""

.. code-block:: LAMMPS

   units lj
   pair_style grav/long 4.0
   pair_coeff * *
   kspace_style ewald/grav 1.0e-8

   pair_style grav/long 4.0 gconst 2.5
   pair_coeff * *

   pair_style hybrid/overlay lj/cut 1.122462 grav/long 5.0
   pair_coeff * * lj/cut 1.0 1.0
   pair_coeff * * grav/long

Used in input scripts:

   .. parsed-literal::

      examples/PACKAGES/ewald_grav/in.jeans
      examples/PACKAGES/ewald_grav/in.fewbody

Description
"""""""""""

.. versionadded:: TBD

The *grav/long* style computes the real-space part of Newtonian gravity
between point masses in a periodic system, for use with the
:doc:`kspace_style ewald/grav <kspace_style>` long-range solver, which
computes the remaining part in reciprocal space.  As for long-range
Coulomb interactions, the 1/r interaction is split into a short-range
part, computed here up to the cutoff, and a smooth long-range part:

.. math::

   E = - \frac{G m_i m_j}{r} \, \textrm{erfc} \left( \alpha r \right) \qquad r < r_c

where :math:`G` is the gravitational constant, :math:`m_i` and
:math:`m_j` are the masses of the two atoms, :math:`r_c` is the cutoff,
and :math:`\alpha` is the Ewald splitting parameter (the G vector)
chosen by the KSpace style, see the :doc:`kspace_modify gewald
<kspace_modify>` command.  The force is the derivative of :math:`E`,

.. math::

   F = - \frac{G m_i m_j}{r^2} \left[ \textrm{erfc} \left( \alpha r \right)
       + \frac{2 \alpha r}{\sqrt{\pi}} \exp \left( - \alpha^2 r^2 \right) \right]

directed along the line between the atoms, i.e. attractive.  The method
follows :ref:`(Hernquist) <Hernquist1991>`, where :math:`E` is the
real-space term of Eq. (2.11).

.. note::

   Eq. (2.14b) of :ref:`(Hernquist) <Hernquist1991>` is missing the factor
   :math:`r = |\mathbf{x} - \mathbf{n} L|` in front of the exponential
   term, likely a misprint.  LAMMPS uses the correct form shown above.

The masses are the per-type masses set with the :doc:`mass <mass>`
command or in the data file, or the per-atom masses of atom styles
that define them, e.g. :doc:`atom_style sphere <atom_style>`.  The
atoms do not need a charge.

The gravitational constant can be set with the *gconst* keyword.  If it
is not given, the value for the current :doc:`units <units>` is used. The
*ewald/grav* style takes G from this pair style, so that both parts of
the interaction always use the same value.

Gravity in the real units of atomic and molecular systems is extremely
weak.  For astrophysical or cosmological applications it is often most
convenient to use :doc:`units lj <units>` with G = 1 and to scale the
results to physical units afterwards.

The interaction has no softening: two point masses that come close
experience the full :math:`1/r^2` force.  To avoid close encounters, a
suitably small timestep can be used, or a short-range repulsion can be
added with :doc:`pair_style hybrid/overlay <pair_hybrid>`, as in the
last example above.

The pair_coeff command takes no coefficients other than the atom types,
and all type pairs use the global cutoff.  The interactions between
atoms that are bonded to each other are weighted by the Coulomb factors
of the :doc:`special_bonds <special_bonds>` command.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Mixing is not relevant for this pair style, since it has no per-type
coefficients.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style writes its information, including the gravitational
constant, to :doc:`binary restart files <restart>`, so pair_style and
pair_coeff commands do not need to be specified in an input script that
reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style requires a gravity KSpace style, currently
:doc:`kspace_style ewald/grav <kspace_style>`.

This pair style is part of the KSPACE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`kspace_modify <kspace_modify>`, :doc:`pair_style coul/long <pair_coul>`

Default
"""""""

gconst = value for the current units, see the table above

----------

.. _Hernquist1991:

**(Hernquist)** L. Hernquist, F. R. Bouchet, and Y. Suto, ApJS, 75, 231 (1991).
